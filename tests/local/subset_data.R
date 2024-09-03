library(dplyr)

# Load full data
raw_data <- data.table::fread(input='../TS_ALUNO_34EM.csv',
                          integer64='character',
                          skip=0,  
                          nrow=-1, 
                          na.strings = "", 
                          showProgress = TRUE)

# Randomly select 320 schools
set.seed(22408)
#set.seed(78244) good for RJ
#set.seed(846274) good for country-wise
saeb <- raw_data |> 
  dplyr::filter(!is.na(INSE_ALUNO), IN_PROFICIENCIA_MT == 1, ID_UF == 33
  ) |>
  dplyr::select(ID_ESCOLA, IN_PUBLICA, INSE_ALUNO, PROFICIENCIA_MT, ID_LOCALIZACAO)

saeb <- saeb |> 
  dplyr::rename(
    public = IN_PUBLICA,
    student_ses = INSE_ALUNO,
    math_proficiency = PROFICIENCIA_MT,
    location = ID_LOCALIZACAO
  )

selected_schools <- sample(unique(saeb$ID_ESCOLA), 320, replace = FALSE)

saeb <- saeb[ID_ESCOLA %in% selected_schools, ]

saeb$school_id <- 0
k <- 0
for( i in unique(saeb$ID_ESCOLA) ) {
  k <- k+1
  saeb[saeb$ID_ESCOLA == i, "school_id"] <- k
}

saeb <- saeb |> 
  dplyr::select(school_id, public, location, student_ses, math_proficiency)

save(saeb, file = "data/saeb.rda")
