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
                   column(12, plotlyOutput("summary_outliers_calc", height = "1110"))
                 ),
                 br(),
                 fluidRow(
                   column(12, plotlyOutput("summary_outliers_phase", height = "500px"))
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
                          downloadButton("download_outlier_calcs", "Download Calculation Outliers"),
                          downloadButton("download_outlier_phase", "Download Phase Outliers"),
                          downloadButton("download_outlier_firstval", "Download First-Value Outliers"),
                          downloadButton("download_outlier_all", "Download All Outliers"),
                          br(), br(),
                          DT::dataTableOutput("outlier_data_table")
                   )
                 )
        )
      )
    )
  )
)


# Server
server <- function(input, output, session) {
  data_result <- reactiveVal(NULL)
  segments_data <- data.frame(
    start_time = c(0, 35.5, 86.62, 154.75, 205.88),
    end_time = c(35.5, 86.62, 154.75, 205.88, 265.42),
    condition = c("Basal 2.8mM Glucose", "16.7mM Glucose", "Oligomycin (5uM)", "FCCP (3uM)", "Rotenone/Antimycin A (5uM)"),
    color = c("Blue", "Green", "Indigo", "Grey", "Magenta"),
    text_position = c(18, 60, 120, 180, 235)
  )
  colors <- c("#7774FF", "#92C797", "#C274FF", "#85A3BE", "#F974FF")

  observeEvent(input$process, {
    req(input$ocr_file, input$dna_file)
    df <- rio::import(input$ocr_file$datapath)
    dna <- rio::import(input$dna_file$datapath)

    result <- process_ocr_data(df, dna, input$num_islets)
    data_result(result)
  })

  observeEvent(data_result(), {
    df <- data_result()$ocr_traces
    updateSelectizeInput(session, "selected_donor",
                         choices = sort(unique(df$donor_id)),
                         server = TRUE)
  })

  observeEvent(input$reset_donor, {
    updateSelectizeInput(session, "selected_donor", selected = "")
  })

  # Processed Data tab
  # Calculations data tab
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

  # Traces data tab
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

  output$outlier_data_table <- DT::renderDataTable({
    req(data_result())

    calc_df <- data_result()$outlier_calcs %>%
      dplyr::filter(FLAG_outliercalc_whennormdna) %>%
      dplyr::mutate(outlier_type = "Calculation")

    phase_df <- data_result()$outlier_phase %>%
      dplyr::filter(FLAG_outlierphase_whennormdnabaseline) %>%
      dplyr::mutate(outlier_type = "Phase")

    firstval_df <- data_result()$outlier_baseline %>%
      dplyr::mutate(
        outlier_type = "First value",
        flagged_firstval = FLAG_firstvallessthat50_raw | FLAG_firstvalmorethan450_raw
      ) %>%
      dplyr::filter(flagged_firstval)

    combined <- dplyr::bind_rows(
      calc_df,
      phase_df,
      firstval_df
    )

    DT::datatable(
      combined,
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

  output$download_outlier_all <- downloadHandler(
    filename = function() paste0("ocr_outlier_all_", Sys.Date(), ".csv"),
    content = function(file) {
      calc_df <- data_result()$outlier_calcs %>%
        dplyr::filter(FLAG_outliercalc_whennormdna) %>%
        dplyr::mutate(outlier_type = "Calculation")

      phase_df <- data_result()$outlier_phase %>%
        dplyr::filter(FLAG_outlierphase_whennormdnabaseline) %>%
        dplyr::mutate(outlier_type = "Phase")

      firstval_df <- data_result()$outlier_baseline %>%
        dplyr::mutate(
          outlier_type = "First value",
          flagged_firstval = FLAG_firstvallessthat50_raw | FLAG_firstvalmorethan450_raw
        ) %>%
        dplyr::filter(flagged_firstval)

      combined <- dplyr::bind_rows(calc_df, phase_df, firstval_df)
      write.csv(combined, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)

