library(rmarkdown)

print("Generating report.")

# first command line argument is the prefix of the source files, second genome length
args <- commandArgs(TRUE)
source = args[1]
genome_length = as.integer(args[2])

# generate from Rmarkdown file html report with processed input data
rmarkdown::render("./report.Rmd", params = list(source_path = source, genome_len = genome_length), output_file = paste(source, "-report.html", sep = ""))
