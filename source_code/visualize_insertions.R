# visualize inserted bases in the text file
# dot means no change, concrete base represents inserted sequence

library("data.table")
print("Visualizing insertions")

# first command line argument is the prefix of the source files, second genome length
args <- commandArgs(TRUE)
source = args[1]
genome_len = as.integer(args[2])

write_to_file <- function(x) {
  position = as.integer(x[2])
  sequence = x[4]
  seq_length = as.integer(x[3])
  line = paste(paste(replicate(position - 1, "."), collapse = ""), sequence, paste(replicate(max(c(0, genome_len - position - seq_length + 1)), "."), collapse = ""), sep="")
  write(line, file = paste(source, "-insertion_visualization.txt", sep=""),
        ncolumns = 1,
        append = TRUE, sep = "\n")
}

table <- fread(paste(source, "-insertions.txt", sep=""), header = TRUE)
if (nrow(table) == 0) {
  file.create(paste(source, "-insertion_visualization.txt", sep = ""))
  stop("There are no insertions in the reads.")
}


insertions <- table[table$End == "I"]
if (nrow(insertions) == 0) {
  file.create(paste(source, "-insertion_visualization.txt", sep = ""))
  print("There are no insertions in the reads.")
} else {

insertions <- insertions[order(insertions$Position, insertions$Length, insertions$Sequence), ]


apply(insertions, 1, write_to_file)
}
