############################################################
# CRC Explorer: Colorectal Cancer Biomarker Analysis Tool
# A Shiny app for analyzing gene expression data in colorectal cancer
# Includes built-in sample data for immediate testing
############################################################

# Load required libraries
library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)
library(survival)
library(survminer)
library(heatmaply)
library(RColorBrewer)

# Function to generate sample data
generate_sample_data <- function() {
  # Set seed for reproducibility
  set.seed(123)
  
  # Number of patients
  n_patients <- 100
  
  # Create patient IDs
  patient_ids <- paste0("P", sprintf("%03d", 1:n_patients))
  
  # Generate expression values for common CRC biomarkers
  gene_exp_data <- data.frame(
    patient_id = patient_ids,
    APC = rnorm(n_patients, mean = 5.2, sd = 1.2),
    TP53 = rnorm(n_patients, mean = 6.8, sd = 1.5),
    KRAS = rnorm(n_patients, mean = 4.5, sd = 1.8),
    MLH1 = rnorm(n_patients, mean = 3.7, sd = 1.1),
    BRAF = rnorm(n_patients, mean = 2.9, sd = 1.3),
    PIK3CA = rnorm(n_patients, mean = 5.6, sd = 1.4)
  )
  
  # Round expression values to 2 decimal places
  gene_exp_data[, 2:7] <- round(gene_exp_data[, 2:7], 2)
  
  # Create stages with different distributions
  stages <- sample(c("Stage I", "Stage II", "Stage III", "Stage IV"), 
                   size = n_patients, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
  
  # Create treatment response categories
  response_categories <- sample(c("Complete Response", "Partial Response", "Stable Disease", "Progressive Disease"),
                                size = n_patients, replace = TRUE, prob = c(0.2, 0.3, 0.3, 0.2))
  
  # Generate survival data
  time_to_event <- round(rexp(n_patients, rate = 1/36) + 1, 1)  # Mean survival of 36 months
  time_to_event <- pmin(time_to_event, 60)  # Cap at 60 months
  
  # Event indicator (1 = event occurred, 0 = censored)
  event <- rbinom(n_patients, size = 1, prob = 0.6)
  
  # Create survival status category
  survival_status <- ifelse(event == 1, "Deceased", "Alive")
  
  # Create clinical data frame
  clinical_data <- data.frame(
    patient_id = patient_ids,
    age = round(rnorm(n_patients, mean = 62, sd = 10)),
    gender = sample(c("Male", "Female"), n_patients, replace = TRUE),
    stage = stages,
    response = response_categories,
    time_to_event = time_to_event,
    event = event,
    survival_status = survival_status
  )
  
  # Adjust TP53 based on stage
  stage_factor <- as.factor(clinical_data$stage)
  stage_effect <- as.numeric(factor(stage_factor, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))) * 0.5
  
  gene_exp_data$TP53 <- gene_exp_data$TP53 + stage_effect - 1
  
  # Adjust KRAS based on survival status
  survival_effect <- ifelse(clinical_data$survival_status == "Deceased", 1.2, 0)
  gene_exp_data$KRAS <- gene_exp_data$KRAS + survival_effect
  
  # Round expression values again after adjustments
  gene_exp_data[, 2:7] <- round(gene_exp_data[, 2:7], 2)
  
  # Return the data frames
  return(list(
    gene_expression = gene_exp_data,
    clinical_data = clinical_data
  ))
}

# Generate sample data
sample_data <- generate_sample_data()

# UI Component
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel(div(
    img(src = "crc_logo.png", height = 60),
    "CRC Explorer: Colorectal Cancer Biomarker Analysis Tool",
    style = "display: flex; align-items: center; gap: 15px;"
  )),
  
  sidebarLayout(
    # Sidebar with input options
    sidebarPanel(
      width = 3,
      h4("Input Options"),
      
      # File uploads
      h5("Upload Files"),
      fileInput("gene_exp_file", "Gene Expression Data (CSV)",
                accept = c("text/csv", "text/comma-separated-values", ".csv")),
      fileInput("clinical_file", "Clinical Data (CSV)",
                accept = c("text/csv", "text/comma-separated-values", ".csv")),
      
      # Sample data option
      checkboxInput("use_sample_data", "Use Sample Data", value = TRUE),
      
      # Biomarker selection
      h5("Biomarker Selection"),
      selectInput("biomarker_select", "Select Biomarkers",
                  choices = c("APC", "TP53", "KRAS", "MLH1", "BRAF", "PIK3CA"),
                  multiple = TRUE,
                  selected = c("APC", "TP53", "KRAS")),
      
      # Analysis options
      h5("Analysis Options"),
      selectInput("group_var", "Group by:",
                  choices = c("Cancer Stage", "Treatment Response", "Survival Status"),
                  selected = "Cancer Stage"),
      
      selectInput("analysis_type", "Analysis:", 
                  choices = c("Expression Comparison", "Survival Analysis", "Correlation Analysis"),
                  selected = "Expression Comparison"),
      
      actionButton("analyze_btn", "Analyze", class = "btn-primary btn-block")
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Box Plot", 
                 plotlyOutput("boxplot"),
                 br(),
                 textOutput("stat_result"),
                 downloadButton("download_boxplot", "Download Results")
        ),
        tabPanel("Survival", 
                 plotOutput("survival_plot", height = "500px"),
                 br(),
                 verbatimTextOutput("survival_stats"),
                 downloadButton("download_survival", "Download Results")
        ),
        tabPanel("Heatmap", 
                 plotlyOutput("heatmap_plot"),
                 br(),
                 downloadButton("download_heatmap", "Download Results")
        ),
        tabPanel("Data Explorer", 
                 DT::dataTableOutput("data_table"),
                 br(),
                 downloadButton("download_data", "Download Data")
        ),
        tabPanel("About", 
                 br(),
                 h3("About CRC Explorer"),
                 p("CRC Explorer is an interactive web-based tool for analyzing gene expression patterns of colorectal cancer biomarkers across different patient cohorts."),
                 p("Built using R Shiny, CRC Explorer enables researchers to visualize biomarker expression patterns, perform survival analyses, and identify correlations between multiple biomarkers in an intuitive interface."),
                 h4("Features:"),
                 tags$ul(
                   tags$li("Intuitive interface requiring minimal training to operate effectively"),
                   tags$li("Compare expression of multiple biomarkers across different patient subgroups"),
                   tags$li("Perform survival analysis to evaluate prognostic potential of biomarkers"),
                   tags$li("Visualize correlations between different biomarkers"),
                   tags$li("Export publication-ready figures and data tables")
                 ),
                 h4("Sample Data:"),
                 p("The tool includes a synthetic dataset of 100 colorectal cancer patients with expression values for six common biomarkers (APC, TP53, KRAS, MLH1, BRAF, PIK3CA) and clinical information."),
                 h4("Citation:"),
                 p("If you use CRC Explorer in your research, please cite:"),
                 p(em("Maleki S. (2024). CRC Explorer: An Interactive Web-Based Tool for Colorectal Cancer Biomarker Analysis. Bilkent University."))
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Reactive data for gene expression
  gene_data <- reactive({
    if (input$use_sample_data) {
      return(sample_data$gene_expression)
    } else {
      req(input$gene_exp_file)
      df <- read.csv(input$gene_exp_file$datapath, header = TRUE)
      return(df)
    }
  })
  
  # Reactive data for clinical information
  clinical_data <- reactive({
    if (input$use_sample_data) {
      return(sample_data$clinical_data)
    } else {
      req(input$clinical_file)
      df <- read.csv(input$clinical_file$datapath, header = TRUE)
      return(df)
    }
  })
  
  # Merged data
  merged_data <- reactive({
    req(gene_data(), clinical_data())
    
    # Assuming both datasets have a common 'patient_id' column
    merged <- merge(gene_data(), clinical_data(), by = "patient_id")
    return(merged)
  })
  
  # Reactive for the selected grouping variable
  grouping_var <- reactive({
    if(input$group_var == "Cancer Stage") {
      return("stage")
    } else if(input$group_var == "Treatment Response") {
      return("response")
    } else {
      return("survival_status")
    }
  })
  
  # Generate box plot 
  output$boxplot <- renderPlotly({
    req(merged_data(), input$biomarker_select, input$analyze_btn)
    
    # Reshape data for plotting
    plot_data <- merged_data() %>%
      select(patient_id, !!sym(grouping_var()), all_of(input$biomarker_select)) %>%
      pivot_longer(cols = all_of(input$biomarker_select), 
                   names_to = "biomarker", 
                   values_to = "expression")
    
    # Create box plot
    p <- ggplot(plot_data, aes(x = !!sym(grouping_var()), y = expression, fill = biomarker)) +
      geom_boxplot(alpha = 0.7) +
      stat_boxplot(geom = 'errorbar', width = 0.3) +
      theme_bw() +
      labs(x = input$group_var, y = "Expression Level", 
           title = paste("Expression of Selected Biomarkers by", input$group_var)) +
      scale_fill_brewer(palette = "Set1") +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
            plot.title = element_text(hjust = 0.5, size = 14))
    
    ggplotly(p)
  })
  
  # Calculate and show statistics
  output$stat_result <- renderText({
    req(merged_data(), input$biomarker_select, input$analyze_btn)
    
    # Perform ANOVA for each biomarker
    results <- lapply(input$biomarker_select, function(biomarker) {
      formula <- as.formula(paste(biomarker, "~", grouping_var()))
      model <- aov(formula, data = merged_data())
      summary_val <- summary(model)
      p_value <- summary_val[[1]]["Pr(>F)"][[1]][1]
      return(c(biomarker = biomarker, p_value = p_value))
    })
    
    # Combine results
    result_text <- "Statistical Results: "
    for (i in seq_along(results)) {
      result_text <- paste0(result_text, results[[i]][1], " (p-value = ", 
                            round(as.numeric(results[[i]][2]), 4), ")", 
                            ifelse(i < length(results), ", ", ""))
    }
    return(result_text)
  })
  
  # Generate survival plot
  output$survival_plot <- renderPlot({
    req(merged_data(), input$biomarker_select, input$analyze_btn, 
        "time_to_event" %in% colnames(merged_data()), "event" %in% colnames(merged_data()))
    
    if (length(input$biomarker_select) == 0) return(NULL)
    
    # Use first selected biomarker for survival analysis
    biomarker <- input$biomarker_select[1]
    
    # Create high/low expression groups
    surv_data <- merged_data()
    surv_data$expr_group <- ifelse(surv_data[[biomarker]] > median(surv_data[[biomarker]], na.rm = TRUE), 
                                   "High Expression", "Low Expression")
    
    # Fit survival model
    fit <- survfit(Surv(time_to_event, event) ~ expr_group, data = surv_data)
    
    # Plot survival curves
    ggsurvplot(fit, data = surv_data, pval = TRUE, risk.table = TRUE,
               xlab = "Time (months)", ylab = "Survival Probability",
               title = paste("Survival Analysis by", biomarker, "Expression"),
               palette = c("#E7B800", "#2E9FDF"),
               ggtheme = theme_bw())
  })
  
  # Calculate and show survival statistics
  output$survival_stats <- renderPrint({
    req(merged_data(), input$biomarker_select, input$analyze_btn,
        "time_to_event" %in% colnames(merged_data()), "event" %in% colnames(merged_data()))
    
    if (length(input$biomarker_select) == 0) return(NULL)
    
    # Use first selected biomarker
    biomarker <- input$biomarker_select[1]
    
    # Create high/low expression groups
    surv_data <- merged_data()
    surv_data$expr_group <- ifelse(surv_data[[biomarker]] > median(surv_data[[biomarker]], na.rm = TRUE), 
                                   "High Expression", "Low Expression")
    
    # Cox model
    cox_model <- coxph(Surv(time_to_event, event) ~ expr_group, data = surv_data)
    summary(cox_model)
  })
  
  # Generate heatmap
  output$heatmap_plot <- renderPlotly({
    req(merged_data(), input$biomarker_select, input$analyze_btn)
    
    if (length(input$biomarker_select) < 2) {
      return(NULL)
    }
    
    # Extract biomarker expression data
    heatmap_data <- merged_data() %>%
      select(all_of(input$biomarker_select))
    
    # Calculate correlation
    cor_matrix <- cor(heatmap_data, use = "pairwise.complete.obs")
    
    # Create heatmap
    heatmaply(cor_matrix, 
              colors = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                          "#F4A582", "#FDDBC7", "#FFFFFF", 
                                          "#D1E5F0", "#92C5DE", "#4393C3", 
                                          "#2166AC", "#053061"))(100),
              main = "Correlation Heatmap of Selected Biomarkers",
              xlab = "Biomarkers", 
              ylab = "Biomarkers",
              showticklabels = c(TRUE, TRUE),
              margins = c(50, 50, 40, 40))
  })
  
  # Show data table
  output$data_table <- renderDT({
    req(merged_data(), input$analyze_btn)
    
    # Show only relevant columns
    display_data <- merged_data() %>%
      select(patient_id, all_of(input$biomarker_select), !!sym(grouping_var()))
    
    datatable(display_data, 
              options = list(pageLength = 10, scrollX = TRUE, 
                             autoWidth = TRUE, searching = TRUE),
              rownames = FALSE)
  })
  
  # Download handlers
  output$download_boxplot <- downloadHandler(
    filename = function() {
      paste("CRC_boxplot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(merged_data(), input$biomarker_select)
      
      # Reshape data
      plot_data <- merged_data() %>%
        select(patient_id, !!sym(grouping_var()), all_of(input$biomarker_select)) %>%
        pivot_longer(cols = all_of(input$biomarker_select), 
                     names_to = "biomarker", 
                     values_to = "expression")
      
      # Create plot
      p <- ggplot(plot_data, aes(x = !!sym(grouping_var()), y = expression, fill = biomarker)) +
        geom_boxplot(alpha = 0.7) +
        stat_boxplot(geom = 'errorbar', width = 0.3) +
        theme_bw() +
        labs(x = input$group_var, y = "Expression Level", 
             title = paste("Expression of Selected Biomarkers by", input$group_var)) +
        scale_fill_brewer(palette = "Set1") +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              plot.title = element_text(hjust = 0.5, size = 14))
      
      ggsave(file, plot = p, width = 10, height = 7, dpi = 300)
    }
  )
  
  output$download_survival <- downloadHandler(
    filename = function() {
      paste("CRC_survival_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(merged_data(), input$biomarker_select,
          "time_to_event" %in% colnames(merged_data()), "event" %in% colnames(merged_data()))
      
      biomarker <- input$biomarker_select[1]
      
      # Create high/low expression groups
      surv_data <- merged_data()
      surv_data$expr_group <- ifelse(surv_data[[biomarker]] > median(surv_data[[biomarker]], na.rm = TRUE), 
                                     "High Expression", "Low Expression")
      
      # Fit model
      fit <- survfit(Surv(time_to_event, event) ~ expr_group, data = surv_data)
      
      # Create plot
      p <- ggsurvplot(fit, data = surv_data, pval = TRUE, risk.table = TRUE,
                      xlab = "Time (months)", ylab = "Survival Probability",
                      title = paste("Survival Analysis by", biomarker, "Expression"),
                      palette = c("#E7B800", "#2E9FDF"),
                      ggtheme = theme_bw())
      
      # Save plot
      ggsave(file, plot = print(p), width = 10, height = 8, dpi = 300)
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste("CRC_heatmap_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(merged_data(), input$biomarker_select)
      
      # Extract data
      heatmap_data <- merged_data() %>%
        select(all_of(input$biomarker_select))
      
      # Calculate correlation
      cor_matrix <- cor(heatmap_data, use = "pairwise.complete.obs")
      
      # Create heatmap
      p <- heatmaply(cor_matrix, 
                     colors = colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                                                 "#F4A582", "#FDDBC7", "#FFFFFF", 
                                                 "#D1E5F0", "#92C5DE", "#4393C3", 
                                                 "#2166AC", "#053061"))(100),
                     file = file,
                     width = 10, height = 8,
                     main = "Correlation Heatmap of Selected Biomarkers",
                     xlab = "Biomarkers", 
                     ylab = "Biomarkers",
                     showticklabels = c(TRUE, TRUE),
                     margins = c(50, 50, 40, 40))
    }
  )
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste("CRC_data_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(merged_data())
      
      # Get displayed data
      display_data <- merged_data() %>%
        select(patient_id, all_of(input$biomarker_select), !!sym(grouping_var()))
      
      # Write to CSV
      write.csv(display_data, file, row.names = FALSE)
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)