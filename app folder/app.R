library(shiny)
library(tidyverse)
library(plotly)
library(purrr)
library(rio)
library(pracma)
library(gtsummary)
library(ggplot2)
library(shinyBS)

# Define function
source("process_ocr.R")

# plotting function
segments_data <- data.frame(
  start_time = c(0, 35.5, 86.62, 154.75, 205.88),
  end_time = c(35.5, 86.62, 154.75, 205.88, 265.42),
  condition = c("2.8mM Glucose", "16.7mM Glucose", "Oligomycin (5uM)",
                "FCCP (3uM)", "Rotenone/Antimycin A (5uM)"),
  color = c("#7774FF", "#92C797", "#C274FF", "#85A3BE", "#F974FF"),
  text_position = c(18, 60, 120, 180, 235)
)
make_interactive_plot <- function(df, norm_type, selected_donor = NULL) {

  # Filter data for the normalisation type
  plot_data <- df %>%
    dplyr::filter(normalisation_type == norm_type) %>%
    dplyr::mutate(donor_id = as.factor(donor_id)) %>%
    dplyr::arrange(donor_id, time)

  # Set y-axis label based on normalisation type
  y_axis_title <- dplyr::case_when(
    norm_type == "Raw" ~ "OCR (pmol/min)",
    norm_type == "Norm DNA" ~ "OCR (pmol/min/µg DNA)",
    norm_type == "Norm DNA and Baseline" ~ "OCR (fold change from baseline)",
    TRUE ~ "OCR"
  )

  #  donor selection/highlighting
  shared_data <- highlight_key(plot_data, key = ~donor_id)

  # condition bands so they stay fixed
  rect_shapes <- lapply(seq_len(nrow(segments_data)), function(i) {
    list(
      type = "rect",
      x0 = segments_data$start_time[i],
      x1 = segments_data$end_time[i],
      xref = "x",
      y0 = 0.92,
      y1 = 1.00,
      yref = "paper",
      fillcolor = segments_data$color[i],
      line = list(width = 0),
      opacity = 0.25,
      layer = "below"
    )
  })

  # Condition labels centered in each band
  segment_annotations <- lapply(seq_len(nrow(segments_data)), function(i) {
    list(
      x = segments_data$text_position[i],
      y = 0.96,
      xref = "x",
      yref = "paper",
      text = segments_data$condition[i],
      showarrow = FALSE,
      xanchor = "center",
      yanchor = "middle",
      font = list(size = 11, color = "black")
    )
  })

  p <- plot_ly(data = shared_data) %>%
    add_lines(
      x = ~time,
      y = ~mean,
      split = ~donor_id,
      line = list(width = 1.5, color = "magenta"),
      hoverinfo = "text",
      text = ~paste0(
        "Donor: ", donor_id,
        "<br>Time: ", round(time, 2),
        "<br>Value: ", round(mean, 3)
      ),
      showlegend = FALSE
    ) %>%
    layout(
      title = list(
        text = norm_type,
        x = 0.5,
        xanchor = "center",
        font = list(size = 16, family = "Arial", color = "black")
      ),
      margin = list(l = 75, r = 40, t = 80, b = 70),
      plot_bgcolor = "rgba(0,0,0,0)",
      paper_bgcolor = "rgba(0,0,0,0)",
      xaxis = list(
        title = "Time (minutes)",
        showgrid = FALSE,
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "black",
        tickcolor = "black",
        ticks = "outside",
        tickfont = list(color = "black", family = "Arial"),
        titlefont = list(color = "black", family = "Arial", size = 12)
      ),
      yaxis = list(
        title = y_axis_title,
        showgrid = FALSE,
        zeroline = FALSE,
        showline = TRUE,
        linecolor = "black",
        tickcolor = "black",
        ticks = "outside",
        tickfont = list(color = "black", family = "Arial"),
        titlefont = list(color = "black", family = "Arial", size = 12),
        type = "linear"
      ),
      shapes = rect_shapes,
      annotations = segment_annotations
    ) %>%
    highlight(
      on = "plotly_click",
      off = "plotly_doubleclick",
      persistent = FALSE,
      opacityDim = 0.15,
      selected = attrs_selected(line = list(color = "black", width = 3)),
      defaultValues = if (!is.null(selected_donor) && nzchar(selected_donor)) {
        list(selected_donor)
      } else {
        NULL
      }
    )

  return(p)
}
make_interactive_outlier_calc_plot <- function(df) {
  nice_names <- c(
    calc_atp_resp = "ATP-linked respiration",
    calc_atp_resp_gluc = "ATP-linked (glucose-stim.)",
    calc_basal_resp = "Basal respiration",
    calc_bhi = "BHI",
    calc_bhi_gluc = "BHI (glucose stim)",
    calc_gluc_stim_oci = "Glucose-stimulated OCI (max/basal)",
    calc_gluc_stim_sparecap = "Spare capacity (glucose-stim.)",
    calc_stim_gluc_resp = "Stimulated glucose respiration",
    calc_max_resp = "Max respiration (FCCP)",
    calc_max_gluc_resp = "Max glucose respiration (16.7mM)",
    calc_nonmito_oc = "Non-mitochondrial OCR",
    calc_proton_leak = "Proton leak",
    calc_spare_cap = "Spare capacity"
  )

  df2 <- df %>%
    dplyr::mutate(
      measurement_label = dplyr::recode(measurement, !!!nice_names, .default = measurement),
      outlier_flag = ifelse(FLAG_outliercalc_whennormdna, "Outlier", "Not outlier")
    )

  panels <- split(df2, df2$measurement_label)
  panel_names <- names(panels)

  panel_plots <- lapply(panel_names, function(m) {
    d <- panels[[m]]
    set.seed(1)
    d$x_jit <- runif(nrow(d), 0.85, 1.15)

    plot_ly(
      data = d,
      x = ~x_jit,
      y = ~value,
      type = "scatter",
      mode = "markers",
      color = ~outlier_flag,
      colors = c("Not outlier" = "forestgreen", "Outlier" = "red3"),
      text = ~paste0(
        "Donor: ", donor_id,
        "<br>Metric: ", measurement_label,
        "<br>Value: ", round(value, 2),
        "<br>Status: ", outlier_flag
      ),
      hoverinfo = "text",
      showlegend = FALSE,
      marker = list(size = 7, opacity = 0.8)
    ) %>%
      layout(
        xaxis = list(
          title = "",
          showticklabels = FALSE,
          zeroline = FALSE,
          range = c(0.75, 1.25)
        ),
        yaxis = list(title = "", zeroline = FALSE),
        margin = list(l = 45, r = 10, t = 25, b = 20)
      )
  })

  subplot(
    panel_plots,
    nrows = ceiling(length(panel_plots) / 3),
    shareX = FALSE,
    shareY = FALSE,
    titleX = FALSE,
    titleY = FALSE,
    margin = 0.04
  ) %>%
    layout(
      title = "Calculation outliers",
      annotations = lapply(seq_along(panel_names), function(i) {
        list(
          text = panel_names[i],
          x = ((i - 1) %% 3 + 0.5) / 3,
          y = 1 - (floor((i - 1) / 3)) / ceiling(length(panel_names) / 3) + 0.02,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          xanchor = "center",
          font = list(size = 11)
        )
      })
    )
}


make_interactive_outlier_phase_plot <- function(df) {
  df2 <- df %>%
    dplyr::mutate(
      outlier_flag = ifelse(FLAG_outlierphase_whennormdnabaseline, "Outlier", "Not outlier")
    )

  panels <- split(df2, df2$phase)
  panel_names <- names(panels)

  panel_plots <- lapply(panel_names, function(ph) {
    d <- panels[[ph]]
    set.seed(1)
    d$x_jit <- runif(nrow(d), 0.85, 1.15)

    plot_ly(
      data = d,
      x = ~x_jit,
      y = ~auc,
      type = "scatter",
      mode = "markers",
      color = ~outlier_flag,
      colors = c("Not outlier" = "forestgreen", "Outlier" = "red3"),
      text = ~paste0(
        "Donor: ", donor_id,
        "<br>Phase: ", phase,
        "<br>AUC: ", round(auc, 2),
        "<br>Status: ", outlier_flag
      ),
      hoverinfo = "text",
      showlegend = FALSE,
      marker = list(size = 7, opacity = 0.8)
    ) %>%
      layout(
        xaxis = list(
          title = "",
          showticklabels = FALSE,
          zeroline = FALSE,
          range = c(0.75, 1.25)
        ),
        yaxis = list(title = "", zeroline = FALSE),
        margin = list(l = 45, r = 10, t = 25, b = 20)
      )
  })

  subplot(
    panel_plots,
    nrows = ceiling(length(panel_plots) / 3),
    shareX = FALSE,
    shareY = FALSE,
    titleX = FALSE,
    titleY = FALSE,
    margin = 0.05
  ) %>%
    layout(
      title = "Phase AUC outliers",
      annotations = lapply(seq_along(panel_names), function(i) {
        list(
          text = panel_names[i],
          x = ((i - 1) %% 3 + 0.5) / 3,
          y = 1 - (floor((i - 1) / 3)) / ceiling(length(panel_names) / 3) + 0.03,
          xref = "paper",
          yref = "paper",
          showarrow = FALSE,
          xanchor = "center",
          font = list(size = 11)
        )
      })
    )
}

# UI
ui <- fluidPage(

  titlePanel("OCR Seahorse Data Analysis"),

  sidebarLayout(
    sidebarPanel(
      fileInput("ocr_file", "Upload Seahorse file", accept = c(".csv", ".xlsx", ".xls")),
      fileInput("dna_file", "Upload DNA file", accept = c(".csv", ".xls", ".xlsx")),
      numericInput("num_islets", "Islets per well", value = 70),
      actionButton("process", "Run Analysis"),
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("OCR traces",
                 # Container with search and reset bar at top
                 fluidRow(
                   column(12,
                          div(style = "padding: 15px; background: rgba(255,255,255,0.9); border-radius: 8px; box-shadow: 0px 0px 10px #ccc;",
                              selectizeInput("selected_donor", "Search Donor ID:",
                                             choices = NULL,
                                             multiple = FALSE,
                                             options = list(placeholder = 'Type a donor ID to highlight...')),
                              actionButton("reset_donor", "Reset Highlight", icon = icon("eraser"))
                          )
                   )
                 ),
                 br(),

                 # Graphs below the search bar
                 fluidRow(
                   column(12, plotlyOutput("plot_raw", height = "400px"))
                 ),
                 fluidRow(
                   column(12, plotlyOutput("plot_norm_dna", height = "400px"))
                 ),
                 fluidRow(
                   column(12, plotlyOutput("plot_norm_dna_baseline", height = "400px"))
                 )
        ),
        tabPanel("Outlier Summary",
                 fluidRow(
                   column(12,
                          tags$b("Triplicate-level QC message"),
                          br(),
                          verbatimTextOutput("triplicate_outlier_message_summary"),
                          uiOutput("triplicate_download_ui")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(12, plotlyOutput("summary_outliers_calc", height = "1100px"))
                 ),
                 br(),
                 fluidRow(
                   column(12, plotlyOutput("summary_outliers_phase", height = "450px"))
                 )
        ),
        tabPanel("Calculations Data",
                 fluidRow(
                   column(12,
                          downloadButton("download_calculations", "Download Calculations Data"),
                          br(), br(),
                          DT::dataTableOutput("calculations_data_table")
                   )
                 )
        ),
        tabPanel("Traces Data",
                 fluidRow(
                   column(12,
                          downloadButton("download_traces", "Download Traces Data"),
                          br(), br(),
                          DT::dataTableOutput("traces_data_table")
                   )
                 )
        ),
        tabPanel("Outlier Data",
                 fluidRow(
                   column(12,
                          tags$b("Triplicate-level QC message"),
                          br(),
                          verbatimTextOutput("triplicate_outlier_message")
                   )
                 ),
                 br(),

                 h4("Calculation Outliers"),
                 fluidRow(
                   column(12,
                          downloadButton("download_outlier_calcs", "Download Calculation Outliers"),
                          br(), br(),
                          DT::dataTableOutput("outlier_calcs_table")
                   )
                 ),
                 br(),

                 h4("Phase AUC Outliers"),
                 fluidRow(
                   column(12,
                          downloadButton("download_outlier_phase", "Download Phase Outliers"),
                          br(), br(),
                          DT::dataTableOutput("outlier_phase_table")
                   )
                 ),
                 br(),

                 h4("First-Value Outliers"),
                 fluidRow(
                   column(12,
                          downloadButton("download_outlier_firstval", "Download First-Value Outliers"),
                          br(), br(),
                          DT::dataTableOutput("outlier_firstval_table")
                   )
                 ),
                 br(),

                 h4("Triplicate Outliers"),
                 fluidRow(
                   column(12,
                          uiOutput("triplicate_download_ui_2"),
                          uiOutput("triplicate_table_ui")
                   )
                 )
        )
      )
    )
  )
)


server <- function(input, output, session) {
  data_result <- reactiveVal(NULL)

  observeEvent(input$process, {
    req(input$ocr_file, input$dna_file)

    df <- rio::import(input$ocr_file$datapath)
    dna <- rio::import(input$dna_file$datapath)

    result <- process_ocr_data(df, dna, input$num_islets)
    data_result(result)
  })

  observeEvent(data_result(), {
    req(data_result())

    df <- data_result()$ocr_traces
    updateSelectizeInput(
      session,
      "selected_donor",
      choices = sort(unique(df$donor_id)),
      server = TRUE
    )
  })

  observeEvent(input$reset_donor, {
    updateSelectizeInput(session, "selected_donor", selected = "")
  })

  # ----------------------------
  # Main data tables
  # ----------------------------
  output$calculations_data_table <- DT::renderDataTable({
    req(data_result())

    df <- data_result()$ocr_calculations

    DT::datatable(
      df,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        scrollY = "600px",
        scrollX = TRUE,
        paging = FALSE,
        ordering = TRUE,
        filter = "top"
      ),
      rownames = FALSE
    )
  })

  output$download_calculations <- downloadHandler(
    filename = function() paste0("ocr_calculations_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(data_result()$ocr_calculations, file, row.names = FALSE)
    }
  )

  output$traces_data_table <- DT::renderDataTable({
    req(data_result())

    df <- data_result()$ocr_traces

    DT::datatable(
      df,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        scrollY = "600px",
        scrollX = TRUE,
        paging = FALSE,
        ordering = TRUE,
        filter = "top"
      ),
      rownames = FALSE
    )
  })

  output$download_traces <- downloadHandler(
    filename = function() paste0("ocr_traces_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(data_result()$ocr_traces, file, row.names = FALSE)
    }
  )

  # ----------------------------
  # Interactive plots
  # ----------------------------
  output$summary_outliers_calc <- renderPlotly({
    req(data_result())
    make_interactive_outlier_calc_plot(data_result()$outlier_calcs)
  })

  output$summary_outliers_phase <- renderPlotly({
    req(data_result())
    make_interactive_outlier_phase_plot(data_result()$outlier_phase)
  })

  output$plot_raw <- renderPlotly({
    req(data_result())
    make_interactive_plot(
      df = data_result()$ocr_traces,
      norm_type = "Raw",
      selected_donor = input$selected_donor
    )
  })

  output$plot_norm_dna <- renderPlotly({
    req(data_result())
    make_interactive_plot(
      df = data_result()$ocr_traces,
      norm_type = "Norm DNA",
      selected_donor = input$selected_donor
    )
  })

  output$plot_norm_dna_baseline <- renderPlotly({
    req(data_result())
    make_interactive_plot(
      df = data_result()$ocr_traces,
      norm_type = "Norm DNA and Baseline",
      selected_donor = input$selected_donor
    )
  })

  # ----------------------------
  # Outlier downloads
  # ----------------------------
  output$download_outlier_calcs <- downloadHandler(
    filename = function() paste0("ocr_outlier_calculations_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- data_result()$outlier_calcs %>%
        dplyr::filter(FLAG_outliercalc_whennormdna)

      write.csv(df, file, row.names = FALSE)
    }
  )

  output$download_outlier_phase <- downloadHandler(
    filename = function() paste0("ocr_outlier_phase_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- data_result()$outlier_phase %>%
        dplyr::filter(FLAG_outlierphase_whennormdnabaseline)

      write.csv(df, file, row.names = FALSE)
    }
  )

  output$download_outlier_firstval <- downloadHandler(
    filename = function() paste0("ocr_outlier_firstvalue_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- data_result()$outlier_baseline %>%
        dplyr::mutate(
          flagged_firstval = FLAG_firstvallessthat50_raw | FLAG_firstvalmorethan450_raw
        ) %>%
        dplyr::filter(flagged_firstval)

      write.csv(df, file, row.names = FALSE)
    }
  )

  output$download_outlier_triplicate <- downloadHandler(
    filename = function() paste0("ocr_outlier_triplicate_", Sys.Date(), ".csv"),
    content = function(file) {
      df <- data_result()$outlier_triplicate %>%
        dplyr::arrange(donor_id, replicate, time)

      write.csv(df, file, row.names = FALSE)
    }
  )

  # ----------------------------
  # Triplicate outlier message
  # ----------------------------
  output$triplicate_outlier_message <- renderText({
    req(data_result())

    trip_df <- data_result()$outlier_triplicate

    if (is.null(trip_df) || nrow(trip_df) == 0) {
      "No outliers were identified within triplicate measurements."
    } else {
      ids <- trip_df %>%
        dplyr::mutate(donor_replicate = paste0(donor_id, "_", replicate)) %>%
        dplyr::pull(donor_replicate) %>%
        unique() %>%
        sort()

      paste0(
        "Potential outliers were identified within triplicate measurements for: ",
        paste(ids, collapse = ", "),
        ". See the table below for the corresponding timepoints. ",
        "Replicates flagged as potential outliers are not removed from the dataset; ",
        "these flags are provided for transparency and quality control."
      )
    }
  })

  output$triplicate_outlier_message_summary <- renderText({
    req(data_result())

    trip_df <- data_result()$outlier_triplicate

    if (is.null(trip_df) || nrow(trip_df) == 0) {
      "No outliers were identified within triplicate measurements."
    } else {
      ids <- trip_df %>%
        dplyr::mutate(donor_replicate = paste0(donor_id, "_", replicate)) %>%
        dplyr::pull(donor_replicate) %>%
        unique() %>%
        sort()

      paste0(
        "Potential outliers were identified within triplicate measurements for: ",
        paste(ids, collapse = ", "),
        ". See the table below for the corresponding timepoints. ",
        "Replicates flagged as potential outliers are not removed from the dataset; ",
        "these flags are provided for transparency and quality control."
      )
    }
  })

  # ----------------------------
  # Conditional triplicate UI
  # ----------------------------
  output$triplicate_download_ui <- renderUI({
    req(data_result())

    trip_df <- data_result()$outlier_triplicate

    if (is.null(trip_df) || nrow(trip_df) == 0) {
      return(NULL)
    }

    downloadButton("download_outlier_triplicate", "Download Triplicate Outliers")
  })

  output$triplicate_download_ui_2 <- renderUI({
    req(data_result())

    trip_df <- data_result()$outlier_triplicate

    if (is.null(trip_df) || nrow(trip_df) == 0) {
      return(NULL)
    }

    tagList(
      downloadButton("download_outlier_triplicate", "Download Triplicate Outliers"),
      br(), br()
    )
  })

  output$triplicate_table_ui <- renderUI({
    req(data_result())

    trip_df <- data_result()$outlier_triplicate

    if (is.null(trip_df) || nrow(trip_df) == 0) {
      return(NULL)
    }

    DT::dataTableOutput("outlier_triplicate_table")
  })

  # ----------------------------
  # Outlier tables
  # ----------------------------
  output$outlier_calcs_table <- DT::renderDataTable({
    req(data_result())

    df <- data_result()$outlier_calcs %>%
      dplyr::filter(FLAG_outliercalc_whennormdna)

    DT::datatable(
      df,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        scrollX = TRUE,
        pageLength = 20
      ),
      rownames = FALSE
    )
  })

  output$outlier_phase_table <- DT::renderDataTable({
    req(data_result())

    df <- data_result()$outlier_phase %>%
      dplyr::filter(FLAG_outlierphase_whennormdnabaseline)

    DT::datatable(
      df,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        scrollX = TRUE,
        pageLength = 20
      ),
      rownames = FALSE
    )
  })

  output$outlier_firstval_table <- DT::renderDataTable({
    req(data_result())

    df <- data_result()$outlier_baseline %>%
      dplyr::mutate(
        flagged = FLAG_firstvallessthat50_raw | FLAG_firstvalmorethan450_raw
      ) %>%
      dplyr::filter(flagged)

    DT::datatable(
      df,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        scrollX = TRUE,
        pageLength = 20
      ),
      rownames = FALSE
    )
  })

  output$outlier_triplicate_table <- DT::renderDataTable({
    req(data_result())

    df <- data_result()$outlier_triplicate
    req(!is.null(df), nrow(df) > 0)

    df <- df %>%
      dplyr::arrange(donor_id, replicate, time)

    DT::datatable(
      df,
      extensions = "Buttons",
      options = list(
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        scrollX = TRUE,
        pageLength = 20
      ),
      rownames = FALSE
    )
  })
}

shinyApp(ui, server)



