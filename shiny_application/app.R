#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(ggplot2)
source('Functions_score.R')

#seq = 'atgcatcatcaccatcaccatGCACCGGTTCGTAGCCTGAATTGTACCCTGCGTGATAGCCAGCAGAAAAGCCTGGTTATGAGCGGTCCGTATGAACTGAAAGCACTGCATCTGCAGGGTCAGGATATGGAACAGCAGGTTGTTTTTAGCATGAGCTTTGTTCAGGGTGAAGAAAGCAATGATAAAATTCCGGTTGCACTGGGTCTGAAAGAAAAGAATCTGTATCTGTCATGTGTTCTGAAAGATGATAAACCGACCCTGCAGCTGGAAAGCGTTGATCCGAAATCGTATCCGAAAAAGAAAATGGAAAAACGTTTTGTGTTTAATAAAATTGAAATTAATAATAAACTGGAATTTGAAAGCGCACAGTTTCCGAATTGGTATATTAGCACCAGCCAGGCAGAAAATATGCCGGTGTTTCTGGGTGGCACCAAAGGTGGTCAGGATATTACCGATTTTACCATGCAGTTTGTTAGCAGCTAA'

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("SCTR"),
    
    titlePanel(h4("SCTR (sequence context for TCG recoding) is an application used to calculate score, design and redesign sequence for UAA incorporation to TCG sense codon.")),
    titlePanel(h5("For Score calculation, a nucleic acid sequence in length of at least 15 is needed. Each codon will be changed to TCG for the score calculation except the first two and the last two codons.")),
    titlePanel(h5("For Sequence design from AA, an amino acid sequence in length of 5 is needed and the mid amino acid will be changed to U for incorporation.")),
    titlePanel(h5("For Sequence redesign from NC, a nucleic acid sequence in length of 15 is needed and the mid codon will be changed to TCG for incorporation.")),
    titlePanel(h5("After calculation, the dataset (as .csv) and plot (as .pdf) can be downloaded when click the buttons.")),
    titlePanel(h5("#the weights file was updated at 2023-6-13.#")),
    
    sidebarLayout(
        sidebarPanel(
          selectInput('calc_type', 'I want to do', choices = c('Score calculation','Sequence design from AA','Sequence redesign from NC')),
          textAreaInput('seq',h3('Sequence'), placeholder = 'Enter sequence', ),
          actionButton('x1','Calculate'),
          downloadButton('downloadData','Download Dataset'),
          downloadButton('downloadPlot','Download Plot')
        ),

        mainPanel(
          plotOutput('score_plot',brush=brushOpts(id='plot_brush')),
          verbatimTextOutput("brush_info"),
          DTOutput('score_result')
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    df_calc = function(ct,seq_i){
      if (ct == 'Score calculation'){
        return(seq_to_score_matrix(seq_i))
      }
      else if (ct == 'Sequence design from AA'){
        return(aaseq_redesign_score_matrix(seq_i))
      }
      else if (ct == 'Sequence redesign from NC'){
        return(ncseq_redesign_score_matrix(seq_i))
      }
    }
    
    dfplot_calc = function(ct){
      if (ct == 'Score calculation'){
        p = ggplot(df(), aes(x=AA_site, y=Score)) + geom_point(size=2) + theme_bw()
        return(p)
      }
      else if (ct == 'Sequence design from AA'){
        p = ggplot(df(), aes(x=Seq, y=Score)) + geom_point(size=2) + xlab('Seq') + theme_bw() +
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank())
        return(p)
      }
      else if (ct == 'Sequence redesign from NC'){
        p = ggplot(df(), aes(x=Seq, y=Score)) + geom_point(size=2) + xlab('Seq') + theme_bw() + 
          theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank())
        return(p)
      }
    }
    
    df <- eventReactive(input$x1, {df_calc(input$calc_type,input$seq)})
    dfplot <- eventReactive(input$x1, {dfplot_calc(input$calc_type)})

    output$score_result <- renderDT({
      df()
    })
    
    output$score_plot <- renderPlot({
      if (length(input$score_result_rows_selected)){
        if (input$calc_type == 'Score calculation'){
          dfplot() + geom_point(data=df()[input$score_result_rows_selected,],mapping=aes(x=AA_site, y=Score),size=4,color='red')
        }
        else if (input$calc_type == 'Sequence design from AA'){
          dfplot() + geom_point(data=df()[input$score_result_rows_selected,],mapping=aes(x=Seq, y=Score),size=4,color='red')
        }
        else if (input$calc_type == 'Sequence redesign from NC'){
          dfplot() + geom_point(data=df()[input$score_result_rows_selected,],mapping=aes(x=Seq, y=Score),size=4,color='red')
        }
      }
      else{
        dfplot()
      }
    })
    
    output$brush_info <- renderPrint({
      brushedPoints(df(),input$plot_brush)
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste('result', ".csv", sep = "")
      },
      content = function(file) {
        write.csv(df(), file, row.names = FALSE)
      }
    )
    
    output$downloadPlot <- downloadHandler(
      filename = function() {
        paste('result', ".pdf", sep = "")
      },
      content = function(file) {
        ggsave(file, dfplot())
      }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)
