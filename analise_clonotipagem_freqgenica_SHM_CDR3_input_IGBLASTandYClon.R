library(data.table)
library(tidyverse)
library(dplyr)
library(openxlsx)
library(Biostrings)
library(stringr)
library(data.table)
 
#Definindo o caminho das minhas amostras, foi definido pela grnade pasta
meus_arquivos <- "/home/yala/Documentos/YALA/5.CLONOTIPAGEM-ANALISE"

#Cria um dataframe contendo os arquivos de INPUT que sairam do YClon de cada amostra
INPUT <- data.frame(report = sort(c(list.files(path = meus_arquivos, 
                                         recursive = TRUE, 
                                         all.files = TRUE, 
                                         pattern = "clonotyped_report.tsv.gz"))), #todos os arquivos com esse padrao entra na coluna de report
                    clonotipado = sort(c(list.files(path = meus_arquivos, 
                                             recursive = TRUE, 
                                             all.files = TRUE, 
                                             pattern = "clonotyped.tsv.gz")))) #todos os arquivos com esse padrao entra na coluna de clonotipado

#Cria um loop onde todos as análises sao feitas das amostras. 
for(y in 1:nrow(INPUT)){
 print(INPUT[y,])

  #Cria uma variavel REPORT com a coluna report do df INPUT
  REPORT <- fread(paste(meus_arquivos,
                        INPUT$report[y],
                        sep="/"),
                  header = FALSE, sep = "\t")
  # REPORT <- fread(paste(meus_arquivos, 
  #                      "PESADO/834-I-H_db-pass_YClon_clonotyped_report.tsv.gz",
  #                       sep="/"),
  #                 header = FALSE, sep = "\t") 
  #Renomeando os nomes das colunas do report
  colnames(REPORT) <- c("sequence_id", "seq_count","most_common_cdr3","clone_id")
  # Remover tudo que aparece depois do "|" na coluna sequence_id 
  REPORT$sequence_id <- sub("\\|.*", "", REPORT$sequence_id)
  
  #Importar arquivo da Clonotipagem e organizar a tabela 
  CLONOTIPAGEM <- fread(paste(meus_arquivos,
                              INPUT$clonotipado[y],
                              sep="/"),
                  header = TRUE, sep = "\t")
  # CLONOTIPAGEM <- fread(paste(meus_arquivos, 
  #                             "PESADO/834-I-H_db-pass_YClon_clonotyped.tsv.gz", 
  #                             sep="/"),
  #                       header = TRUE, sep = "\t")
  CLONOTIPAGEM <- CLONOTIPAGEM %>% select(sequence_id, v_call, v_identity, j_call, j_identity, cdr3, cdr3_aa, clone_id, clone_seq_count)
  # Ordenar o data frame em ordem decrescente pela coluna 
  CLONOTIPAGEM <- arrange(CLONOTIPAGEM, desc(clone_seq_count)) #nesse momento nao e necessario
  # Remover tudo que aparece depois do "|" na coluna sequence_id no merged_data
  CLONOTIPAGEM$sequence_id <- sub("\\|.*", "", CLONOTIPAGEM$sequence_id)
  #Manter na coluna apenas o 1o anotação e deletar tudo depois da vírgula
  CLONOTIPAGEM$v_call <- c(sub(",(.*)", "",CLONOTIPAGEM$v_call)) #exclui tudo depois da virgula, estuda os alelos
  #CLONOTIPAGEM$v_call <- c(sub("-(.*)", "",CLONOTIPAGEM$v_call)) #exclui tudo depois do -, só estuda familia genica
  CLONOTIPAGEM$j_call <- c(sub(",(.*)", "",CLONOTIPAGEM$j_call)) #exclui tudo depois da virgula, estuda os alelos
  #CLONOTIPAGEM$j_call <- sub("\\*.*", "", CLONOTIPAGEM$j_call) #exclui tudo depois do *, só estuda familia genica
  
  # Junta os df atraves da coluna "sequence_id" que está presente em ambas
  #REPORT <- left_join(REPORT[, names(REPORT) != "clone_id"], CLONOTIPAGEM, by = "sequence_id") #não está funcionado
  REPORT <- merge(REPORT, CLONOTIPAGEM, by = "sequence_id", all.x = TRUE) #CORRIGIR, PQ TEM A MESMA COLUNA 2x

  #Daqui para baixo, é apenas análise dos dados
  
  ##Análise da frequência genica
  #V
  ValorVabsoluto <- table(REPORT$v_call)
  ValorVrelativo <- prop.table(table(REPORT$v_call))
  TabelaGeneV <- data.frame (ValorVabsoluto, ValorVrelativo)
  #J
  ValorJabsoluto <- table(REPORT$j_call)
  ValorJrelativo <- prop.table(table(REPORT$j_call))
  TabelaGeneJ <- data.frame (ValorJabsoluto, ValorJrelativo)
  
  #Hipermutacao somatica
  #Colocar as sequencias de forma rankeada pela coluna do clone_seq_count
  #V
  hiperV <- REPORT %>%
    mutate(v_identity = 100 - v_identity) %>%
    group_by(v_call) %>%
    summarise(hipermutacao_v = mean(v_identity)) 
  #J
  hiperJ <- REPORT %>%
    mutate(j_identity = 100 - j_identity) %>%
    group_by(j_call) %>%
    summarise(hipermutacao_j = mean(j_identity)) 
  
  #Análise do CDR3
  #Tamanho
  REPORT$cdr3_tam <- nchar(REPORT$cdr3_aa)
  #Cria uma tabela de frequencia dos tamanhos do cdr3_aa
  CDR3absoluto <- table(REPORT$cdr3_tam)
  CDR3relativo <- prop.table(table(REPORT$cdr3_tam))
  TAM_CDR3 <- data.frame (CDR3absoluto, CDR3relativo)
  # Criar o histograma da coluna "cdr3_tam"
  #hist(REPORT$cdr3_tam, main = "Distribuição de Frequência do tamanho do CDR3", xlab = "Tamanho CDR3", ylab = "Frequência")
  
  #Análise da composição do CDR3
  # Cria novo data frame apenas com a coluna cdr3_aa da REPORT
  Analise_CDR3 <- data.frame(cdr3_aa = REPORT$cdr3_aa)
  #cria 20 colunas no data frame, e conta quantos caracteres tem na sequencia de aa na coluna cdr3_aa
  amino_acidos <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  #Conta a quantidade de determinados caracteres/aminoacidos que estao nas sequencias do CDR3
  for (aa in amino_acidos) {
    Analise_CDR3[[aa]] <- str_count(REPORT$cdr3_aa, aa)
  }
   COMP_CDR3 <- list()
    for (x in amino_acidos) {
    frequencia <- sum(Analise_CDR3[x]) / sum(REPORT$cdr3_tam)
    COMP_CDR3[[x]] <- frequencia
  }
  #Cria um df com as informações da freq dos aminoacidos
  COMP_CDR3 <- data.frame(amino_acido = names(COMP_CDR3), frequencia = unlist(COMP_CDR3))
  
    # Exportando! Criando arquivos finais
  # Uma unica tabela do excel com 6 abas
   workbook <- createWorkbook()
   addWorksheet(workbook, sheetName = "Freq_V")
   writeData(workbook, sheet = 1, x = TabelaGeneV)
   addWorksheet(workbook, sheetName = "Freq_J")
   writeData(workbook, sheet = 2, x = TabelaGeneJ)
   addWorksheet(workbook, sheetName = "Hiper_V")
   writeData(workbook, sheet = 3, x = hiperV)
   addWorksheet(workbook, sheetName = "Hiper_J")
   writeData(workbook, sheet = 4, x = hiperJ)
   addWorksheet(workbook, sheetName = "Tamanho_CDR3")
   writeData(workbook, sheet = 5, x = TAM_CDR3)
   addWorksheet(workbook, sheetName = "Composicao_CDR3")
   writeData(workbook, sheet = 6, x = COMP_CDR3)

   saveWorkbook(workbook, paste(meus_arquivos,
                                str_replace(INPUT$report[y],
                                                          "report.tsv.gz",
                                                          "analise_resultado.xls"),
                                sep = "/"))

} #fechando loop inicial
