#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

lower <- 3
upper <- 20
# Use formatted labels in facetting
pair_names <- c(
  'random' = "Random pairs",
  'sequence' = "Sequence pairs")
       
# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("Visualizing RMSD distributions in protein structure"),
  fluidRow(
    wellPanel(
    sliderInput(inputId = "flen",
                label = "Select length of fragments",
                value = 5, min = lower, max = upper, step = 1),
    helpText("Swipe to select length of n-mers, from 3 to 20"),
    textOutput('matches')
    #textOutput('test')
    )), #end wellPanel
  
  fluidRow(
    column(8, offset = 0,
           plotOutput('hist')),
    column(3, offset = 0,
           plotOutput('boxplot'))
    
  ) # end Plots row
)


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library("ggplot2")
  library("dplyr")
  library("data.table")
  library("cowplot")
  theme_set(theme_bw(base_size = 20))
  
  
my_data <- reactive({
     
     flen <- input$flen
     out_folder = "out/"
     
     sequence_fn <- paste(out_folder, flen, "-mers_rmsd.txt", sep = "")
     random_fn <- paste(out_folder, flen, "-mers_random.txt", sep = "")
     
     
     print("Reading data")
     sequence <- read.csv(file = sequence_fn, header = T, stringsAsFactors = F)
     sequence$pair <- "sequence"
     random <- read.csv(file = random_fn, header = T, stringsAsFactors = F)
     random$pair <- "random"
     
     
 
     # Data wrangling
     # Mark duplicates
     # Duplicate = two or more entries where f1 and f2.seq_id
     # are identical or
     # entries where f1 and f2.seq_id are the same
     sequence$id <- 1:nrow(sequence)
     sequence$duplicated <- FALSE
     sequence[sequence$f1.seq_id == sequence$f2.seq_id, 12] <- TRUE
     DT <- as.data.table(sequence)
     DT[DT[,.SD[2,.N], by=.(f1.seq_id, f2.seq_id)]$id, 12] <- TRUE
     sequence <- as.data.frame(DT)
     # mark the same number of random entries
     random$duplicated <- c(rep(T, sum(sequence$duplicated)),
                            rep(F, nrow(random) - sum(sequence$duplicated)))
     
     
     # Merge sequence and random data
     full_join(sequence, as.data.table(random))
})
  
ndDT <- reactive({
  
  my_data() %>% filter(duplicated == F)
  
})


     output$hist <- renderPlot({
       
     ggplot(data = filter(my_data(), duplicated == F), aes(rmsd, fill = pair)) +
       geom_histogram(bins = 60) +
       facet_wrap( ~ pair, labeller = as_labeller(pair_names)) +
       guides(fill = FALSE) +
       labs(x = "RMSD", y = "Count")
       
   })
     
     
     output$boxplot <- renderPlot({
       
       ggplot(data = ndDT(), aes(x = pair, y = rmsd, fill = pair)) +
         geom_boxplot() +
         labs(x = "", y = "RMSD") +
         guides(fill = FALSE) +
         scale_x_discrete(labels = pair_names)
       
     })
     
     output$matches <- renderText({
     n <- sum(ndDT()$pair == "sequence")
     paste("Amound of ", input$flen, "-mers detected: ", n, sep = "")
     })
     
     # p_value <- reactive({
     #   wilcox.test(ndDT() %>% filter(pair == "random") %>% select(rmsd) %>% .$rmsd,
     #               ndDT() %>% filter(pair == "sequence") %>% select(rmsd) %>% .$rmsd)$p.value
     # })
     # 
     # output$test <- renderText({
     #   paste("Wilcoxon test's with alternative hypothesis: true location shift is not equal to 0 \n P-value: ", p_value(), sep = "")
     # })
     # 
})

# Run the application 
shinyApp(ui = ui, server = server)

