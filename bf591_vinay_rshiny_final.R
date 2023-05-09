## Author: Vinay Kumar Duggineni
## BU BF591
## Final Project - RShiny application that features exploration of Huntington's disease dataset
library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
#install.packages('ggbeeswarm')
library('ggbeeswarm')
library(colourpicker)
library(DESeq2)
library(genefilter)
library(dplyr)
library(gplots)
library('RColorBrewer')
#install.packages("ggfortify")
library("ggfortify")
#install.packages("shinythemes")
library("shinythemes")
library(readr)
library(stringr)

# To increase the file size capabilities of the file inputs
options(shiny.maxRequestSize=30*1024^2)

# Define UI for application
ui <- fluidPage(
  titlePanel("BF 591 Final Project: Exploration of Huntington's disease dataset"),
  tabsetPanel(
    type = "tabs",
    tabPanel("Sample Information Exploration",
             p("The sample tab allows user to explore the sample information on uploading metadata file"),
             sidebarLayout(
               sidebarPanel(
                 fileInput("Sample_file", "Choose a Sample CSV File", accept = ".csv", placeholder = "Sample Information File"),
                 submitButton("Submit"), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Summary", p("The summary tab provides the statstics summary of sample information.") , tableOutput("summary_table")),
                   tabPanel("Table", p("The table tab provides the sortable sample information table."), DT::dataTableOutput("sample_table")),
                   tabPanel("Plots", p("The plots tab provides a histogram plot for selected characteristic of samples."),
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("hist_x_axis", "X-Axis",
                                             choices = c("PMI" = "Sample_pmi","Age of Death" = "Sample_Age_of_Death", "RIN"="Sample_rin","Sequence Reads"="Sample_mrna_seq")),
                                colourInput("samp_color_1", "Outline color","darkred"),
                                colourInput("samp_color_2", "Fill color color","darkblue"),
                                submitButton("Submit")
                              ),
                              mainPanel(plotOutput("sample_plot"))
                            )
                   )
                 )
               )
             )
    ),
    tabPanel("Counts Matrix Exploration",
             p("The counts tab allows the user to explore the data on uploading a normalized counts matrix file"),
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = "Normalized_counts_matrix", "Choose a Normalized Counts Matrix CSV File", accept = ".csv", placeholder = 'Normalized Counts Matrix File'),
                 sliderInput(inputId = "var_slider", "percentile of variance", min = 0, max = 100, value = 40),
                 sliderInput(inputId = "zero_slider", "no. of samples that are non-zero", min = 0, max = 100, value = 40),
                 submitButton("Submit"), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Table or Text",
                            p("The table tab provides the statistics based on the input parameters, percentile of variance, and number of samples that are non-zero."),
                            tableOutput("counts_table")),
                   tabPanel("Scatter Plot",
                            hr("Here is your scatter plot tailored to your variance input: "),
                            plotOutput("variance_plot"),
                            hr("Here is your scatter plot tailored to your non-zero input: "),
                            plotOutput("zero_plot")),
                   tabPanel("Heat map",
                            p("The heatmap tab provides a heat map of the normalized counts input matrix"), 
                            plotOutput("heatmap")),
                   tabPanel("pca", selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                            selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                            p("The PCA tab provides a pca plot of PC1 and PC2 components"),
                            plotOutput("counts_pca"))
                 )
               )
             )),
    tabPanel("Differential Expression",
             p("The differential expression tab allows the user to explore the results of differential expression on uploading a differential expression analysis results."),
             sidebarLayout(
               sidebarPanel(
                 fileInput("DE_Anaylsis", "Choose a Differential Expression Results CSV File", accept = ".csv",placeholder = 'Differential Expression Data File'),
                 submitButton(text = 'Submit'), width = 3
               ),
               mainPanel(
                 tabsetPanel(
                   type = "tabs",
                   tabPanel("Table",
                            p("The table tab provides sortable table of the differential expression analysis results."),
                            DT::dataTableOutput("table",width = "80%")),
                   tabPanel("Volcano Plot", 
                            p("The volcano tab provides a volcano plot based on the input paramters, namely, X-axis and Y-axis"),
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("x_axis", "X-Axis",
                                             choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                                             selected = "log2FoldChange"),
                                radioButtons("y_axis", "Y-Axis",
                                             choices = c("log2FoldChange" = "log2FoldChange","baseMean" = "baseMean", "lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"),
                                             selected = "padj"),
                                sliderInput(inputId = "padjusted", min = -35, max = 0,label = "Select the magnitude of the p adjusted coloring:", value = -15, step = 1),
                                colourInput("de_base", "Base point color","darkred"),
                                colourInput("de_highlight", "Highlight point color","darkblue"),
                                submitButton("Submit"),width=3
                              ),
                              mainPanel(
                                plotOutput("volcano"))))
                 )
               )
             )
    ),
    tabPanel("Visualization of Individual Gene Expression",
             p("The visualization tab allows the user to visualize data about each gene individually."),
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = 'IGE_countsID', label = 'Choose a Normalized Counts Matrix CSV File',accept = ".csv", placeholder = "Normalized Counts Matrix File"),
                 fileInput(inputId = 'IGE_metaID', label = 'Choose a Sample Information CSV File',accept = ".csv", placeholder = "Sample Metadata File"),
                 selectInput("metachoice", choices = c("Diagnosis"="Sample_diagnosis", "PMI"="Sample_pmi","Age of Death"="Sample_Age_of_Death", "RIN"="Sample_rin", "Sequence Reads"="Sample_mrna_seq"),
                             label = "Select Metadata Category Field", selected = "Diagnosis"),
                 #insert gene search box here
                 textInput("gene", label = "Enter Ensembl ID of gene (For ex: ENSG00000000003.10)"),
                 radioButtons("type_plot", "Choose the type of plot",
                              choices = c("Bar Plot" = "Bar Plot","Scatter Plot" = "Scatter Plot", "Violin Plot" = "Violin Plot" , "Beeswarm Plot" = "Beeswarm Plot")),
                 submitButton(text='Plot', icon = icon("bar-chart-o"), width = '100%')
               ),
               mainPanel(
                 plotOutput("IGE_plot")
               )
             ))
  )
  
)

#  Define server logic 
server <- function(input, output,session) {
  # ----------------------------------------------------------------------------------------------------------------------------
  
  # Tab - 1: Sample Information Exploration
  
  # reading the sample data
  #' @details The purpose of this function is to load in the sample information
  #' .csv file. The user has the ability to input their file 
  sample_data <- reactive({
    req(input$Sample_file)
    data <- read.csv(input$Sample_file$datapath, row.names = 1)
    return(data)
  })
  
  # Summary table 
  #' Here we want to create a dataframe summarizing the sample information file
  #' that the user will input.
  #' @param data This is the data frame that is loaded
  #' @details Column_name This is the column where each variable is mentioned
  #' @details summarized_df This is the dataframe that includes the column 
  #' name, classification,mean and the standard deviation of the variable 
  
  
  summary <- function(data){
    Column_Name <- c("Sample_geo_accession", "Sample_status", "Sample_submission_date", 
                     "Sample_last_update_date", "Sample_type", 
                     "Sample_Channel_Count", "Sample_Organ", "Sample_Organism", 
                     "Sample_Tissue_Type", "Sample_Neurological Diagnosis", 
                     "Sample_pmi", "Sample_Age_of_Death", 'Sample_rin', 
                     "Sample_mrna_seq", "Sample_age_of_onset")
    
    Classification <- c("character", "character", "character", "character", 
                        "character", "double", "character", "character", 
                        "character", "character", "double", "double",
                        "double", "double", "double")
    
    Sample_Channel_Count_average <- mean(data$Sample_Channel_Count)
    Sample_pmi_average <- mean(data$Sample_pmi, na.rm = TRUE)
    Sample_Age_of_Death_average <- mean(data$Sample_Age_of_Death)
    Sample_rin_average <- mean(data$Sample_rin)
    Sample_mrna_seq_average <- mean(data$Sample_mrna_seq)
    Sample_Channel_Count_sd <- sd(data$Sample_Channel_Count)
    Sample_pmi_sd <- sd(data$Sample_pmi, na.rm = TRUE)
    Sample_Age_of_Death_sd <- sd(data$Sample_Age_of_Death)
    Sample_rin_sd <- sd(data$Sample_rin)
    Sample_mrna_seq_sd <- sd(data$Sample_mrna_seq)
    Sample_Neurological_Diagnosis_uniq <- toString(unique(data$Sample_diagnosis))
    
    
    Mean <- c('N/A', 'N/A', 'N/A', 
              'N/A', 'N/A', 
              Sample_Channel_Count_average, 'N/A', 'N/A', 
              'N/A', Sample_Neurological_Diagnosis_uniq, 
              Sample_pmi_average, Sample_Age_of_Death_average, Sample_rin_average, 
              Sample_mrna_seq_average, '41.85')
    
    Standard_deviation <- c('N/A', 'N/A', 'N/A', 
                            'N/A', 'N/A', 
                            Sample_Channel_Count_sd, 'N/A', 'N/A', 
                            'N/A', Sample_Neurological_Diagnosis_uniq, 
                            Sample_pmi_sd, Sample_Age_of_Death_sd, Sample_rin_sd, 
                            Sample_mrna_seq_sd, '10.80')
    
    summarized_df <- data_frame(Column_Name, Classification, Mean, Standard_deviation)
    return(summarized_df)
  }
  
  # plot 
  #' @details  Here the goal is to develop a way for the user to generate a 
  #' histogram plot based on the sample information. 
  #' @param data The df with all of the data that the user has inputted
  #' @param x_name This is x-axis. The user has the ability to choose what 
  #' they would like the x-axis of the plot to be. 
  #' @param color_1 This is the outer color that the user chooses. 
  #' @param color_2 This is the inner color that the user chooses.
  
  plot <- function(data, x_name,color_1,color_2){
    data_histogram <- 
      ggplot (data, aes(x = !!sym(x_name))) +
      geom_histogram(color = color_1, fill = color_2) +
      scale_color_manual(values = c(color_1, color_2))
    theme_bw() +
      theme(legend.position = 'bottom') 
    
    return(data_histogram)
  }
  
  # Summary Table output
  output$summary_table <- renderTable({
    data <- sample_data()
    return(summary(data))
  })
  
  # Sample Data Table output
  output$sample_table <- DT::renderDataTable({ 
    data <- sample_data()
    if (is.null(data)) {
      return()
    }
    isolate({
      return(data)
    })
  })
  
  # plot output
  observeEvent(list(input$hist_x_axis,input$samp_color_1,input$samp_color_2),{
    output$sample_plot <- renderPlot({
      data <- sample_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        return(plot(data,input$hist_x_axis,input$samp_color_1,input$samp_color_2))
      })
    }) 
  })
  
  
  # ----------------------------------------------------------------------------------------------------------------------------
  
  # Tab - 2: Counts Matrix Expression
  
  # reading the counts matrix data
  #' Load in the Data 
  #' @details Here all we want to do is simply load in the counts .csv file 
  #' that the user will input via the "Submit" button. 
  counts_data <- reactive({
    req(input$Normalized_counts_matrix)
    data <- read.csv(input$Normalized_counts_matrix$datapath)
    colnames(data)[1] <- "gene"
    return(data)
  })
  
  #' @details The purpose of this function is to filter the data based 
  #' on what value the of number of samples that are greater the variance
  #' threshold and the zero slider threshhold that the user inputs via the slider 
  #' @param data The counts matrix that the user has imported
  #' @param var_slider The slider that allows the user to choose 
  #' variance percentage 
  #' @param zero_slider The slider that the user uses to include genes with at 
  #' least X samples that are non-zero
  
  
  filtering <- function(data,var_slider,zero_slider){
    counts_matrix <- data[,-1]
    variances <- apply(counts_matrix, 1, var, na.rm = TRUE)
    var_threshold <- quantile(variances, var_slider / 100)
    non_zero_counts <- rowSums(counts_matrix > 0, na.rm = TRUE)
    
    filtered <- counts_matrix[variances >= var_threshold & non_zero_counts >= zero_slider, ]
    
    total_samples <- ncol(counts_matrix)
    total_genes <- nrow(counts_matrix)
    passing_genes <- nrow(filtered)
    not_passing_genes <- total_genes - passing_genes
    table <- data.frame( 
      Information = c("Total Number of Samples","Total Number of Genes","Number of Passing Filter", "Number of Not Passing Filter", "Percent of Passing Filter", "Percent of Not Passing Filter"),
      Statstics = c(as.integer(total_samples),as.integer(total_genes),as.integer(passing_genes),as.integer(not_passing_genes), (passing_genes / total_genes) * 100, (not_passing_genes / total_genes) * 100)
    )
    return(table)
  }
  
  # scatter_plot_variance
  #' @details Here we will build a scatter plot based upon what the user 
  #' inputs via the slider. For plotting the y axis we have taken negative y log
  #' @param counts_matrix the counts matrix that the user inputs 
  #' @param slider The values that the user will input 
  scatter_plot_variance <- function(counts_matrix, slider) {
    county <- counts_matrix[,-1] 
    counts_matrix$variance <- rowVars(as.matrix(county))
    counts_matrix <- counts_matrix[rev(order(counts_matrix$variance)),]
    input_value <- slider/100 
    
    counts_matrix$median <- apply(counts_matrix[,-1], 1, median)
    
    Filter_Method <- (floor(nrow(counts_matrix)*input_value)) > counts_matrix$variance
    
    eruption <- ggplot(data = counts_matrix, aes(x = median, y = -log10(variance)))+
      geom_point(aes(color = Filter_Method)) + 
      theme_bw() + 
      scale_color_manual(values = c('darkred', 'darkblue')) + 
      theme(legend.position = "bottom") 
    
    return(eruption)
  }
  
  # scatter_plot_zero
  #' @details Here we will build a scatter plot based upon what the user 
  #' inputs via the slider. For plotting the y axis we have taken negative y log
  #' @param counts_matrix the counts matrix that the user inputs 
  #' @param slider_zero The values that the user will input 
  scatter_plot_zero <- function(counts_matrix, slider_zero) {
    
    counts_matrix$frequency <- rowSums(counts_matrix != 0)
    counts_matrix$median <- apply(counts_matrix[,-1], 1, median)
    
    Filter_Method <- counts_matrix$frequency > slider_zero
    messi <- ggplot(data = counts_matrix, 
                    aes(x = median, y = -log10(frequency)))+
      geom_point(aes(color = Filter_Method )) + 
      theme_bw() + 
      scale_color_manual(values = c('red', 'darkblue')) + 
      theme(legend.position = "bottom") 
    
    
    return(messi)
    
  }
  
  # heatmap
  #'@details Here we want to produce a heat map post filtering of the counts
  #'matrix that the user has inputted
  #'@param counts_tib This is the counts matrix that the user will input
  #'@param perc_var This is the input from the variance slider 
  #'@param nz_genes This is the input from the zero genes slider 
  #'@param num_colors 
  #'@param palette 
  plot_heatmap <- function(counts_tib, perc_var, nz_genes){
    if (!is.null(input$Normalized_counts_matrix)){
      counts_tib <- as_tibble(counts_tib) %>% mutate(across(starts_with(c("C_", "H_")) & where(is.numeric), ~na_if(.x, 0)))
      counts_tib$no_zeros <- rowSums(is.na(counts_tib))  #make new col, with counts.
      counts_tib <- filter(counts_tib, no_zeros <= nz_genes)
      counts_tib <- log10(counts_tib[,!colnames(counts_tib) %in% c("gene", "no_zeros")]) #exclude the gene names column and log scale the values  
      #produce plot_tib
      plot_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib, MARGIN = 1, FUN = var)) #compute variance to filter the data
      perc_val <- quantile(plot_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
      plot_tib <- filter(plot_tib, variance >= perc_val) #filter the tibble
      hmap <- heatmap.2(as.matrix(plot_tib[-ncol(plot_tib)]), scale = "row", col = brewer.pal(9, "YlOrRd"))
      return(hmap)}
    else{return(NULL)}
  }
  
  # pca
  #'@details Here we will construct a PC1 vs PC2 plot based on the input sequence 
  #'that the user will utilize 
  #'@param data The counts matrix that the user will input 
  #pca <- function(data){
  #  data <- data[, -c(1)]
  #  pca_plotting <- prcomp(data, scale. = TRUE)
  #  plot_project <- autoplot(pca_plotting, data = data, colour = "darkblue")
  #  return(plot_project)
 # }
  pca <- function(counts_tib, var_slider, comp1, comp2){
    if (!is.null(input$Normalized_counts_matrix)){
      #make plot tib-
      filt_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var), .after = gene) #calculate variance for filtering
      perc_val <- quantile(filt_tib$variance, probs = var_slider/100, na.rm = TRUE)   #calculate percentile
      filt_tib <- filter(filt_tib, variance >= perc_val) #filter the tibble
      final_tib <- filt_tib[,-c(1,2)]
      pca_res <- prcomp(t(final_tib), center = FALSE, .scale = TRUE) #transpose the data and perform PCA
      pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
      x <- round(pca_var[as.integer(str_sub(comp1, 3))]*100, 2)
      y <- round(pca_var[as.integer(str_sub(comp2, 3))]*100, 2)
      #produce PCA plot
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca <- ggplot(plot_tib, aes(PC1, PC2))+
        geom_point()+
        labs(title="Princple Component Analysis Plot")+
        xlab(str_c(comp1, x, "% Variance", sep=" "))+
        ylab(str_c(comp2, y, "% Variance", sep=" "))+
        theme_bw()
      return(pca)}
    else{return(NULL)}
  }
  
  
  #info table output
  observeEvent(list(input$var_slider,input$zero_slider),{
    output$counts_table <-  renderTable({
      data <- counts_data()
      return(filtering(data,input$var_slider,input$zero_slider))
    })
  })
  
  
  # variance plot output
  observeEvent(input$var_slider,{
    output$variance_plot <- renderPlot( {
      data <- counts_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        scatter_plot_variance(data, input$var_slider)
      })
    }) 
  })
  
  # zero plot output
  observeEvent(input$zero_slider,{
    output$zero_plot <- renderPlot({
      data <- counts_data()
      if (is.null(data)) {
        return()
      }
      isolate({
        scatter_plot_zero(data, input$zero_slider)
      })
    }) 
  })
  
  #heatmap output
  observeEvent(input$zero_slider,{
    output$heatmap <- renderPlot({
      req(input$Normalized_counts_matrix)
      c_matrix <- counts_data()
      hot <- plot_heatmap(c_matrix, input$var_slider, input$zero_slider) 
      return(hot)
      
    },height = 600, width = 1000)
  })
  
  # pca output
  observeEvent(list(input$var_slider, input$comp1, input$comp2),{
  output$counts_pca <- renderPlot({
    data <- counts_data()
    if (is.null(data)) {
      return()
    }
    isolate({
      pca(data, input$var_slider, input$comp1, input$comp2)
    })
  })
  })
  
  # ----------------------------------------------------------------------------------------------------------------------------
  
  # Tab - 3: Differential Expression
  
  # reading the differential expression data 
  #' Load in the DE Data
  #' 
  #' @details Here all we want to do is simply load in the DE .csv file that the user will 
  #' input via the "Submit" button.
  DE_data <- reactive({
    req(input$DE_Anaylsis)
    data <- read.csv(input$DE_Anaylsis$datapath)
    colnames(data)[1] <- "Gene"
    return(data)
  })
  
  #' Volcano plot
  #'
  #' @param dataf The loaded data frame.
  #' @param x_name The column name to plot on the x-axis
  #' @param y_name The column name to plot on the y-axis
  #' @param de_base One of the colors for the points.
  #' @param de_highlight The other colors for the points. Hexadecimal strings: "#CDC4B5"
  #'
  #' @return A ggplot object of a volcano plot
  volcano_plot <-
    function(dataf, x_name, y_name, slider, color1, color2) {
      if(is.null(dataf))     
        return(NULL)
      dataf<-na.omit(dataf)
      out<- ggplot(dataf,aes(x= !!sym(x_name),y= -log10(!!sym(y_name)),color= !!sym(y_name)<(10^slider)))+ 
        theme_bw()+
        theme(legend.position="bottom")+
        ggtitle('Volcano plot')+
        scale_color_manual(name= paste0("padj < 1 x 10^",toString(slider)), values=c(color1, color2))+
        geom_point()
      return(out)
    }
  
  # Using renderPlot to display the plot
  observeEvent(list(input$x_axis, input$y_axis, input$padjusted, input$de_base, input$de_highlight), {
    output$volcano <- renderPlot( {
      data <- DE_data()
      #data <- data %>% mutate(!!paste("neg_log10", input$y_axis, sep = "_") := -log10(!!sym(input$y_axis)))
      if (is.null(data)) {
        return()
      }
      isolate({
        volcano_plot(data, input$x_axis, input$y_axis,input$padjusted, input$de_base, input$de_highlight)
      })
    }) 
  })
  
  # Using renderDataTable to display the table
  output$table <- DT::renderDataTable({ 
    data <- DE_data()
    if (is.null(data)) {
      return()
    }
    isolate({
      return(data)
    })
  })
  
  
  # ----------------------------------------------------------------------------------------------------------------------------
  
  # Tab - 4:
  load_IGE_counts <- reactive({
    if (!is.null(input$IGE_countsID)){
      counts <- read_csv(input$IGE_countsID$datapath)
      colnames(counts)[1] <- "gene"
      return(counts)}
    else{return(NULL)}
  })
  
  IGE_meta <- reactive({
    if (!is.null(input$IGE_metaID)){
      meta <- read_csv(input$IGE_metaID$datapath)
      return(meta)}
    else{return(NULL)}
  })
  
  # Create the plot
  plot_distro <- function(counts_tib, meta_tib, meta_cat, selectgene,type_plot){
    if (!is.null(input$IGE_metaID) & !is.null(input$IGE_countsID)){
      counts_tib <- column_to_rownames(counts_tib, var = "gene")
      gene_counts <- as.numeric(as.vector(counts_tib[selectgene,]))
      plot_tib <- tibble(Gene_Counts = gene_counts, meta_value = pull(meta_tib, meta_cat))
      if (meta_cat == "Sample_diagnosis"){
        plot <- ggplot(plot_tib, aes(meta_value))+
          geom_bar()+
          theme_bw()+
          labs(title = "Plot of gene counts vs sample diagnosis")+ 
          xlab(meta_cat)
        return(plot)
      }
      else {
        if(type_plot == "Bar Plot"){
          plot <- ggplot(plot_tib, aes(meta_value,Gene_Counts)) +
            geom_bar(stat = "identity") +
            labs(title = paste("Gene Expression of", input$gene, "by", input$category),
                 x = input$category, y = "Normalized Count")+ 
            xlab(meta_cat)
        }
        else if(type_plot == "Scatter Plot"){
          plot <- ggplot(plot_tib, aes(meta_value,Gene_Counts))+
            geom_point()+
            theme_bw()+
            labs(title = str_c("Plot of gene counts vs ", meta_cat))+ 
            xlab(meta_cat)
        }
        else if(type_plot == "Violin Plot"){
          plot <- ggplot(plot_tib, aes(meta_value,Gene_Counts))+
            geom_violin() + geom_jitter(width = 0.1)+
            theme_bw()+
            labs(title = str_c("Plot of gene counts vs ", meta_cat))+ 
            xlab(meta_cat)
        }
        else if(type_plot == "Beeswarm Plot"){
          data_range <- max(gene_counts) - min(gene_counts)
          num_bins <- 30
          bin_width <- data_range / num_bins
          
          plot <- ggplot(plot_tib, aes(meta_value,Gene_Counts))+
            geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5,binwidth = bin_width)+
            theme_bw()+
            labs(title = str_c("Plot of gene counts vs ", meta_cat))+ 
            xlab(meta_cat)
        }
        return(plot)
      }}
    else{return(NULL)}
  }
  observeEvent(list(load_IGE_counts(), IGE_meta(), input$metachoice, input$gene, input$type_plot),{
    output$IGE_plot <- renderPlot({
      plot_distro(load_IGE_counts(), IGE_meta(), input$metachoice, input$gene, input$type_plot)
    }) 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)