library(pdftools)

# ============================================================
# Folder containing PDFs
# ============================================================

pdf_dir <- "dissemination/paper/references/intro_discussion"

# ============================================================
# Custom prefix for text files
# ============================================================

txt_prefix <- "ref_intro_discussion_"

# ============================================================
# Get all PDF files
# ============================================================

pdf_files <- list.files(
  path = pdf_dir,
  pattern = "\\.pdf$",
  full.names = TRUE
)

# ============================================================
# Convert PDFs to TXT
# ============================================================

for (pdf_file in pdf_files) {
  
  # Extract text from all pages
  paper_text <- pdf_text(pdf_file)
  
  # Combine pages into one character string
  paper_text <- paste(
    paper_text,
    collapse = "\n\n"
  )
  
  # Create output filename
  txt_file <- file.path(
    pdf_dir,
    paste0(
      txt_prefix,
      tools::file_path_sans_ext(
        basename(pdf_file)
      ),
      ".txt"
    )
  )
  
  # Write text file
  writeLines(
    text = paper_text,
    con = txt_file
  )
  
  message(
    "Converted: ",
    basename(pdf_file)
  )
}

# ============================================================
# Finished
# ============================================================

message(
  "All PDFs converted to text files in: ",
  pdf_dir
)