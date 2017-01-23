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

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("RMSD of n-mer structures in top100H"),
  wellPanel(
    sliderInput(inputId = "flen",
                label = "Select length of fragments",
                value = 5, min = lower, max = upper, step = 1)
  ), #end wellPanel
  plotOutput('distPlot')
  
)


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  library("ggplot2")
  library("dplyr")
  theme_set(theme_bw(base_size = 20))
  
   output$distPlot <- renderPlot({
     
     flen <- input$flen
     out_folder = "out/"
     #plot_folder = "plots/"
     
     signal_fn <- paste(out_folder, flen, "-mers_rmsd.txt", sep = "")
     random_fn <- paste(out_folder, flen, "-mers_random.txt", sep = "")
     #fragments_fn <- paste(out_folder, flen, "-mer_fragments.csv", sep = "")
     
     print("Reading data")
     signal <- read.csv(file = signal_fn, header = T)
     signal$pair <- "signal"
     random <- read.csv(file = random_fn, header = T)
     random$pair <- "random"
     #fragments <- read.csv(file = fragments_fn, header = T)
     
     
     
     check_intraprotein <- function(my_row) {
       if (my_row[2] == my_row[6]) {
         if(my_row[3] == my_row[7] & my_row[4] == my_row[8]) {
           return(NA)
         } else {
           return(TRUE)
         }
       } else {
         return(FALSE) 
       }
     }
     
     signal$intraprotein <- apply(signal, 1, check_intraprotein)
     #Working with interprotein only
     signal <- filter(signal, intraprotein == FALSE) # keep only interprotein
     
     # Should duplicates be excluded?
     # signal <- signal %>%
     #   group_by(f1.seq_id, f2.seq_id) %>%
     #   filter(row_number() == 1)
     
     idx <- sample(x = 1:nrow(random), size = nrow(signal))
     random <- random[idx, ]
     
     my_data <- full_join(signal, random)
     
     
     print("Plotting histograms")
     ggplot(data = my_data, aes(rmsd, fill = pair)) +
       geom_histogram(bins = 60) +
       facet_wrap( ~ pair)
     #fn <- paste(plot_folder, "duplicates_histogram.png", sep = "")
     
   })
})

# Run the application 
shinyApp(ui = ui, server = server)

