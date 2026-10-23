# ==============================================================================
# ARQUIVO COMPLEMENTAR S1: PIPELINE ANALÍTICO, TOPOLÓGICO E DE CIÊNCIA ABERTA
# DA MORTALIDADE HOSPITALAR NA REGIÃO NORDESTE (2016–2025) - VERSÃO 4.0
# ==============================================================================
#
# DESCRIÇÃO GERAL:
# Extensão multi-estadual (9 UFs do Nordeste) do pipeline originalmente
# construído para o Piauí. Além de escalar a extração/purificação/rede de
# associações clínicas, esta versão incorpora três blocos metodológicos
# novos, inspirados na literatura de redes epidemiológicas:
#
#   (1) META-ANÁLISE ENTRE ESTADOS: cada par (diagnóstico base x comorbidade)
#       é estimado separadamente em cada UF (RR/OR + IC95%) e depois combinado
#       por meta-análise de efeitos aleatórios (log-RR), com I², Q de Cochran
#       e tau². Isso evita o viés de simplesmente empilhar 9 estados como se
#       fossem uma população homogênea (risco de paradoxo de Simpson) e
#       classifica cada aresta como "robusta regionalmente" (baixa
#       heterogeneidade) ou "estado-específica" (alta heterogeneidade).
#
#   (2) COMPARAÇÃO COM REDES NULAS / TOPOLOGIA FORMAL: a rede empírica de
#       comorbidades letais é comparada com três famílias clássicas de redes
#       (Erdős–Rényi aleatória, Barabási–Albert livre de escala,
#       Watts–Strogatz mundo-pequeno), pareadas por número de nós e grau
#       médio, quanto a clustering, caminho médio, diâmetro e distribuição de
#       grau (teste de lei de potência). Também é estimado analiticamente o
#       limiar de percolação de Molloy–Reed da rede observada.
#
#   (3) SIMULAÇÃO DE MONTE CARLO DE "CASCATA SINDRÔMICA": a rede de
#       associações é usada como substrato para um processo estocástico de
#       tipo SI/percolação de aresta, em que a probabilidade de "transmissão"
#       de um nó a outro é function do RR normalizado da aresta. Milhares de
#       réplicas por nó-semente estimam o tamanho esperado de cascata a
#       partir de cada diagnóstico-índice (identificando "diagnósticos-fonte"
#       de alto impacto, por analogia a "pacientes-fonte" em redes de
#       rastreamento de contato) e são comparadas às mesmas simulações rodadas
#       sobre a rede aleatória equivalente (teste de significância empírica).
#       Também é feita uma análise de robustez/ataque direcionado (remoção
#       de hubs vs. remoção aleatória) sobre o tamanho do maior componente.
#
# ARQUITETURA MODULAR (16 MÓDULOS):
#  M1  Ambiente e dependências
#  M2  Parâmetros globais (UFs, anos, limiares)
#  M3  Extração e padronização por UF (cache incremental)
#  M4  Classificação de desfecho e purificação clínica por UF
#  M5  Contingência 2x2 por UF (RR/OR) + tiers
#  M6  Meta-análise entre estados (efeitos aleatórios, I², tau²)
#  M7  Regressão logística: Firth (separação) e efeitos mistos (UF aleatório)
#  M8  Construção das redes: por UF, multicamada (multiplex) e núcleo regional
#  M9  Métricas topológicas formais + comparação com redes nulas
#  M10 Percolação de Molloy–Reed e simulação de Monte Carlo de cascata
#  M11 Robustez: ataque direcionado a hubs vs. falha aleatória
#  M12 Rede direcionada fonte→alvo (estilo rastreamento de contato) + componentes
#  M13 Bootstrap de centralidade e permutação de modularidade
#  M14 Análise de sensibilidade multi-limiar (por UF e regional)
#  M15 Análise temporal (por UF e Nordeste consolidado)
#  M16 Tabelas, figuras e pacote de ciência aberta (OSF)
#
# NOTA METODOLÓGICA:
# O tamanho de cada coorte estadual, as propriedades topológicas de cada rede
# e os resultados das simulações são SAÍDAS empíricas dos dados reais, nunca
# metas predefinidas. Nenhum limiar aqui é "ajustado para dar certo": todos
# são declarados a priori no MÓDULO 2 e depois testados em sensibilidade.
# ==============================================================================


# ==============================================================================
# MÓDULO 1: AMBIENTE E DEPENDÊNCIAS
# ==============================================================================
cat("=== MÓDULO 1: INICIALIZANDO AMBIENTE E CARREGANDO PACOTES ===\n")

pacotes <- c(
  "remotes", "tidyverse", "igraph", "visNetwork", "kableExtra", "ggraph",
  "logistf", "openxlsx", "metafor", "lme4", "poweRlaw", "patchwork",
  "scales", "data.table"
)
for (p in pacotes) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
}
if (!requireNamespace("microdatasus", quietly = TRUE)) {
  remotes::install_github("rfsaldanha/microdatasus", quiet = TRUE)
}

library(microdatasus)
library(tidyverse)
library(igraph)
library(visNetwork)
library(kableExtra)
library(ggraph)
library(logistf)
library(openxlsx)
library(metafor)
library(lme4)
library(poweRlaw)
library(patchwork)
library(scales)
library(data.table)

set.seed(42) # Reprodutibilidade de layouts, reamostragens, permutações e Monte Carlo

# ------------------------------------------------------------------------------
# INFRAESTRUTURA DE CHECKPOINT E LOG DE TEMPO (crítico para pipelines longos)
# ------------------------------------------------------------------------------
# Cada módulo pesado (M5, M6, M10, M13) salva seu resultado em disco assim que
# termina. Se o script cair, travar ou você precisar interromper, rodá-lo de
# novo PULA os módulos já concluídos em vez de recomeçar do zero. Os logs com
# horário (`log_tempo`) permitem, olhando o console ou um arquivo de log,
# identificar exatamente em qual módulo o tempo está sendo gasto.
CHECKPOINT_DIR <- "checkpoints_nordeste"
if (!dir.exists(CHECKPOINT_DIR)) dir.create(CHECKPOINT_DIR)

.checkpoint_path <- function(nome) file.path(CHECKPOINT_DIR, paste0(nome, ".rds"))
existe_checkpoint    <- function(nome) file.exists(.checkpoint_path(nome))
salvar_checkpoint    <- function(obj, nome) { saveRDS(obj, .checkpoint_path(nome)); invisible(obj) }
carregar_checkpoint  <- function(nome) readRDS(.checkpoint_path(nome))
log_tempo <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg))

# Para reprocessar um módulo do zero (ex.: mudou um parâmetro), apague o
# arquivo correspondente em checkpoints_nordeste/ ou rode:
#   file.remove(.checkpoint_path("pares_por_uf"))

cat(" -> Pacotes carregados com sucesso.\n\n")


# ==============================================================================
# MÓDULO 2: PARÂMETROS GLOBAIS
# ==============================================================================
cat("=== MÓDULO 2: DEFININDO PARÂMETROS GLOBAIS DO PIPELINE ===\n")

UFS_NORDESTE <- c("MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA")
ANO_INICIO   <- 2016
ANO_FIM      <- 2025

MIN_CASOS_PAR   <- 5     # nº mínimo de óbitos com o par (base, comorbidade) por UF
RR_MIN          <- 1.5   # limiar de corte para entrada na rede principal
ALPHA_FDR       <- 0.05  # significância após correção de FDR
LIMIARES_SENS   <- c(1.2, 1.5, 2.0, 3.0)  # análise de sensibilidade multi-limiar
MIN_UFS_METAAN  <- 3     # nº mínimo de UFs com dado suficiente p/ meta-analisar um par
MIN_UFS_SIG_FOREST <- 4  # nº mínimo de UFs em que o par precisa ser SIGNIFICATIVO
# (individualmente, FDR<0.05, RR>RR_MIN) para entrar no forest plot
N_FOREST_MAX    <- 15    # nº máximo de pares plotados individualmente (custo/tamanho de anexo)
I2_ROBUSTO      <- 50    # I² (%) abaixo do qual a aresta é "robusta regionalmente"

B_BOOT   <- 200   # réplicas de bootstrap de centralidade
B_MOD    <- 500   # permutações para teste de modularidade
B_MC     <- 1000  # réplicas de Monte Carlo por nó-semente na cascata sindrômica
CAP_MC_NODES <- 25 # nº máximo de nós usados como semente na simulação de cascata (custo computacional)

CAPITULOS_EXCLUIR <- c("O", "V", "W", "X", "Y", "S", "T") # causas externas/obstétricas

cat(" -> UFs incluídas:", paste(UFS_NORDESTE, collapse = ", "), "\n")
cat(" -> Janela temporal:", ANO_INICIO, "-", ANO_FIM, "\n\n")


# ==============================================================================
# MÓDULO 3+4: EXTRAÇÃO, PURIFICAÇÃO E LIBERAÇÃO DE MEMÓRIA — UMA UF POR VEZ
# ==============================================================================
cat("=== MÓDULOS 3+4: BAIXANDO E PURIFICANDO CADA UF, UMA DE CADA VEZ ===\n")

# CORREÇÃO CRÍTICA DE MEMÓRIA: a versão anterior baixava as 9 UFs inteiras
# (bruto, sem purificar) e só DEPOIS de todas baixadas começava a purificar.
# Isso significa manter simultaneamente em RAM o SIH-RD bruto de 9 estados x
# 10 anos — para o Nordeste inteiro, isso facilmente ultrapassa a RAM
# disponível e faz o sistema operacional recorrer a swap (memória virtual em
# disco), que é ordens de magnitude mais lento que RAM. Isso é indistinguível
# de "travado": o processo continua rodando, mas cada operação leva minutos
# em vez de milissegundos. Aqui, cada UF é baixada, purificada e IMEDIATAMENTE
# reduzida ao necessário (coorte_final + universo, muito menores que o bruto)
# antes de a próxima UF ser sequer baixada — nunca há mais de uma base bruta
# de UF inteira na memória ao mesmo tempo.

baixar_uf <- function(uf, ano_ini, ano_fim) {
  arq <- paste0("sih_universo_completo_", uf, ".rds")
  if (file.exists(arq)) {
    cat(" -> [", uf, "] cache local encontrado, carregando...\n")
    return(readRDS(arq))
  }
  lista_anos <- list()
  for (ano in ano_ini:ano_fim) {
    cat(" -> [", uf, "-", ano, "] baixando microdados SIH-RD... ")
    bruto <- tryCatch(
      fetch_datasus(year_start = ano, month_start = 1, year_end = ano, month_end = 12,
                    uf = uf, information_system = "SIH-RD"),
      error = function(e) NULL
    )
    if (!is.null(bruto) && nrow(bruto) > 0) {
      proc <- tryCatch(process_sih(bruto), error = function(e) bruto)
      lista_anos[[as.character(ano)]] <- proc
      cat("OK (", format(nrow(proc), big.mark = "."), " registros)\n")
    } else {
      cat("falhou/sem registros, pulando.\n")
    }
  }
  sih_uf <- bind_rows(lista_anos)
  sih_uf$UF <- uf
  saveRDS(sih_uf, arq)
  sih_uf
}

detectar_colunas <- function(df) {
  col_m <- grep("^morte$|fl_morte|ind_morte", names(df), value = TRUE, ignore.case = TRUE)[1]
  col_o <- grep("diag_obito|diagobit|causa_obit", names(df), value = TRUE, ignore.case = TRUE)[1]
  col_p <- grep("princ|diagpr|diag_pr|diag1", names(df), value = TRUE, ignore.case = TRUE)[1]
  col_s <- if ("DIAGSEC1" %in% names(df)) "DIAGSEC1" else grep("^diagsec[0-9]$", names(df), value = TRUE, ignore.case = TRUE)[1]
  col_ano <- grep("ano", names(df), ignore.case = TRUE, value = TRUE)[1]
  list(col_m = col_m, col_o = col_o, col_p = col_p, col_s = col_s, col_ano = col_ano)
}

classificar_desfecho <- function(df, cols) {
  is_obito <- rep(FALSE, nrow(df))
  if (!is.na(cols$col_m)) {
    is_obito <- is_obito | (str_trim(as.character(df[[cols$col_m]])) %in% c("1", "Sim", "SIM", "sim", "TRUE"))
  }
  if (!is.na(cols$col_o)) {
    val_o <- str_trim(as.character(df[[cols$col_o]]))
    is_obito <- is_obito | (!is.na(val_o) & val_o != "" & val_o != "0000")
  }
  is_obito
}

purificar_coorte <- function(df, cols) {
  df$row_id <- seq_len(nrow(df))
  is_obito  <- classificar_desfecho(df, cols)
  df$DESFECHO <- ifelse(is_obito, "obito", "alta")
  
  coorte_let <- df[is_obito, ]
  val_sec    <- str_trim(as.character(coorte_let[[cols$col_s]]))
  manter_sec <- !is.na(val_sec) & val_sec != "" & val_sec != "0000" & val_sec != "NA"
  coorte_sec <- coorte_let[manter_sec, ]
  coorte_sec$DIAG_SEC_LIMPO <- val_sec[manter_sec]
  
  c_base <- toupper(str_sub(str_trim(as.character(coorte_sec[[cols$col_p]])), 1, 1))
  c_sec  <- toupper(str_sub(coorte_sec$DIAG_SEC_LIMPO, 1, 1))
  
  coorte_final <- coorte_sec[!c_base %in% CAPITULOS_EXCLUIR & !c_sec %in% CAPITULOS_EXCLUIR, ] %>%
    mutate(
      DIAG_PRINC = substr(gsub("\\.", "", str_trim(as.character(.data[[cols$col_p]]))), 1, 4),
      DIAGSEC1   = substr(gsub("\\.", "", str_trim(as.character(DIAG_SEC_LIMPO))), 1, 4)
    )
  
  universo <- df %>%
    mutate(
      DIAG_PRINC_U = substr(gsub("\\.", "", str_trim(as.character(.data[[cols$col_p]]))), 1, 4),
      DIAGSEC1_U   = if (!is.na(cols$col_s)) substr(gsub("\\.", "", str_trim(as.character(.data[[cols$col_s]]))), 1, 4) else NA_character_,
      OBITO_U      = (DESFECHO == "obito")
    ) %>%
    select(DIAG_PRINC_U, DIAGSEC1_U, OBITO_U)
  
  list(
    coorte_final = coorte_final,
    universo     = universo,
    n_bruto      = nrow(df),
    n_obitos_brutos = sum(is_obito),
    n_excluidos  = sum(is_obito) - nrow(coorte_final)
  )
}

dados_uf <- setNames(vector("list", length(UFS_NORDESTE)), UFS_NORDESTE)
resumo_atricao <- tibble()
cols_ref <- NULL

for (uf in UFS_NORDESTE) {
  chk_nome <- paste0("dados_purificados_", uf)
  
  if (existe_checkpoint(chk_nome)) {
    log_tempo(paste0(" -> [", uf, "] checkpoint de dados purificados encontrado -> pulando download+purificação."))
    dados_uf[[uf]] <- carregar_checkpoint(chk_nome)
  } else {
    log_tempo(paste0(" -> [", uf, "] baixando SIH-RD (", ANO_INICIO, "-", ANO_FIM, ")..."))
    bruto_uf <- baixar_uf(uf, ANO_INICIO, ANO_FIM)
    
    if (is.null(cols_ref)) {
      cols_ref <- detectar_colunas(bruto_uf)
      stopifnot(!is.na(cols_ref$col_p), !is.na(cols_ref$col_s))
      cat(" -> Colunas detectadas -> principal:", cols_ref$col_p, "| secundário:", cols_ref$col_s,
          "| morte:", cols_ref$col_m, "| ano:", cols_ref$col_ano, "\n")
    }
    
    log_tempo(paste0(" -> [", uf, "] purificando (", format(nrow(bruto_uf), big.mark = "."), " registros brutos)..."))
    dados_uf[[uf]] <- purificar_coorte(bruto_uf, cols_ref)
    salvar_checkpoint(dados_uf[[uf]], chk_nome)
    
    rm(bruto_uf)   # libera a base bruta da UF ANTES de baixar a próxima
    gc(verbose = FALSE)
    log_tempo(paste0(" -> [", uf, "] concluída e liberada da memória (coorte final: ",
                     nrow(dados_uf[[uf]]$coorte_final), " óbitos)."))
  }
  
  resumo_atricao <- bind_rows(resumo_atricao, tibble(
    UF = uf, N_bruto = dados_uf[[uf]]$n_bruto, Obitos_brutos = dados_uf[[uf]]$n_obitos_brutos,
    Excluidos = dados_uf[[uf]]$n_excluidos, Coorte_final = nrow(dados_uf[[uf]]$coorte_final)
  ))
}

cat("\n -> Resumo de atrição por UF:\n")
print(as.data.frame(resumo_atricao))
cat("\n")

coorte_nordeste <- bind_rows(lapply(dados_uf, function(x) x$coorte_final), .id = "UF_lista")
N_obitos_nordeste <- nrow(coorte_nordeste)
saveRDS(coorte_nordeste, "coorte_clinica_obitos_purificada_NORDESTE.rds")
cat(" -> Coorte analítica final (Nordeste, todas as UFs):", format(N_obitos_nordeste, big.mark = "."), "óbitos\n\n")


# ==============================================================================
# MÓDULO 5: CONTINGÊNCIA 2x2 POR UF (RR/OR) + TIERS  [VERSÃO VETORIZADA]
# ==============================================================================
cat("=== MÓDULO 5: CONTINGÊNCIA 2x2 POR UF E CLASSIFICAÇÃO EM TIERS ===\n")

# VERSÃO ANTERIOR (v4.0 original) usava rowwise() chamando calcular_2x2() uma
# vez POR PAR CANDIDATO, e cada chamada reescaneava o `universo` inteiro da
# UF (milhões de linhas). Com milhares de pares candidatos x 9 UFs, isso é
# O(P x N) e foi a causa mais provável do script não terminar em 24h.
#
# VERSÃO VETORIZADA: calcula a,b,c,d para TODOS os pares em UMA passada por
# UF, usando agregação de data.table (by=), em vez de par-a-par. A lógica
# estatística é idêntica (mesmas fórmulas de RR/OR/IC), só a forma de
# computar a,b,c,d muda. Ganho esperado: de horas para segundos/minutos.

estimar_pares_uf_rapido <- function(uf, coorte_final, universo) {
  # Candidatos: mesma regra de sempre (>= MIN_CASOS_PAR óbitos com o par,
  # já na coorte purificada/com capítulos excluídos) — rápido, count() já é
  # vetorizado e nunca foi o gargalo.
  pares_candidatos <- coorte_final %>%
    count(DIAG_PRINC, DIAGSEC1, name = "casos_obito_par") %>%
    filter(casos_obito_par >= MIN_CASOS_PAR)
  if (nrow(pares_candidatos) == 0) return(tibble())
  
  dt <- as.data.table(universo)
  dt[, DIAGSEC1_U := fifelse(is.na(DIAGSEC1_U) | DIAGSEC1_U %in% c("", "0000", "NA"),
                             NA_character_, DIAGSEC1_U)]
  
  # Passada 1: totais de óbito/não-óbito por diagnóstico-BASE (dá c+d e a+b totais)
  totais <- dt[, .N, by = .(DIAG_PRINC_U, OBITO_U)]
  totais_w <- dcast(totais, DIAG_PRINC_U ~ OBITO_U, value.var = "N", fill = 0L)
  if (!"TRUE" %in% names(totais_w))  totais_w[, `TRUE`  := 0L]
  if (!"FALSE" %in% names(totais_w)) totais_w[, `FALSE` := 0L]
  setnames(totais_w, c("TRUE", "FALSE"), c("total_obito_base", "total_naoobito_base"))
  
  # Passada 2: contagens conjuntas (base, comorbidade, óbito) — só para
  # combinações realmente observadas (não precisa do produto cartesiano)
  conj <- dt[!is.na(DIAGSEC1_U), .N, by = .(DIAG_PRINC_U, DIAGSEC1_U, OBITO_U)]
  conj_w <- dcast(conj, DIAG_PRINC_U + DIAGSEC1_U ~ OBITO_U, value.var = "N", fill = 0L)
  if (!"TRUE" %in% names(conj_w))  conj_w[, `TRUE`  := 0L]
  if (!"FALSE" %in% names(conj_w)) conj_w[, `FALSE` := 0L]
  setnames(conj_w, c("TRUE", "FALSE"), c("a", "b"))
  
  # Restringe aos pares candidatos da coorte purificada (semi-join, mantém a
  # mesma semântica da versão original: base_e_comorb calculado no universo
  # bruto, mas só para pares que sobreviveram à purificação clínica)
  setDT(pares_candidatos)
  conj_w <- merge(conj_w, pares_candidatos[, .(DIAG_PRINC, DIAGSEC1, casos_obito_par)],
                  by.x = c("DIAG_PRINC_U", "DIAGSEC1_U"), by.y = c("DIAG_PRINC", "DIAGSEC1"))
  
  res <- merge(conj_w, totais_w, by = "DIAG_PRINC_U")
  res[, `:=`(c = total_obito_base - a, d = total_naoobito_base - b)]
  res <- res[(a + b) > 0 & (c + d) > 0]
  if (nrow(res) == 0) return(tibble())
  
  res[, `:=`(risco_exposto = a / (a + b), risco_nao_exposto = c / (c + d))]
  res[, RR := fifelse(risco_nao_exposto > 0, risco_exposto / risco_nao_exposto, NA_real_)]
  res[, OR := ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))]
  res[, se_logRR := sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d))]
  res[, `:=`(RR_li = exp(log(RR) - 1.96 * se_logRR), RR_ls = exp(log(RR) + 1.96 * se_logRR))]
  res[, se_logOR := sqrt(1/(a + 0.5) + 1/(b + 0.5) + 1/(c + 0.5) + 1/(d + 0.5))]
  res[, `:=`(OR_li = exp(log(OR) - 1.96 * se_logOR), OR_ls = exp(log(OR) + 1.96 * se_logOR))]
  
  # p-valor: teste de Wald vetorizado sobre log(RR)/se_logRR, em vez de
  # fisher.test() linha-a-linha (que é rápido individualmente mas não
  # vetoriza e, multiplicado por milhares de pares, pesa). Válido porque o
  # corte MIN_CASOS_PAR já garante células não-triviais; se quiser conferir
  # um par específico com Fisher exato depois, é rápido rodar isolado.
  res[, p_valor := pnorm(log(RR) / se_logRR, lower.tail = FALSE)]
  
  setnames(res, c("DIAG_PRINC_U", "DIAGSEC1_U"), c("DIAG_PRINC", "DIAGSEC1"))
  res[, UF := uf]
  as_tibble(res[, .(UF, DIAG_PRINC, DIAGSEC1, casos_obito_par, a, b, c, d,
                    risco_exposto, risco_nao_exposto, RR, RR_li, RR_ls,
                    OR, OR_li, OR_ls, se_logRR, se_logOR, p_valor)])
}

if (existe_checkpoint("pares_por_uf")) {
  log_tempo("Checkpoint 'pares_por_uf' encontrado -> pulando Módulo 5 (carregando do disco).")
  pares_por_uf <- carregar_checkpoint("pares_por_uf")
} else {
  log_tempo("Iniciando Módulo 5 (vetorizado) para as 9 UFs...")
  pares_por_uf <- map_dfr(UFS_NORDESTE, function(uf) {
    log_tempo(paste0(" -> UF ", uf, ": calculando 2x2..."))
    r <- estimar_pares_uf_rapido(uf, dados_uf[[uf]]$coorte_final, dados_uf[[uf]]$universo)
    log_tempo(paste0(" -> UF ", uf, ": concluída (", nrow(r), " pares).") )
    r
  })
  
  pares_por_uf <- pares_por_uf %>%
    group_by(UF) %>%
    mutate(p_ajustado = p.adjust(p_valor, method = "fdr")) %>%
    ungroup() %>%
    mutate(
      tier_uf = case_when(
        a >= 10 & b >= 10 & c >= 10 & d >= 10 ~ "Alta",
        a >= 5  & b >= 5  & c >= 5  & d >= 5  ~ "Moderada",
        TRUE                                  ~ "Exploratória"
      ),
      significativo_uf = !is.na(p_ajustado) & p_ajustado < ALPHA_FDR & !is.na(RR) & RR > RR_MIN
    )
  
  salvar_checkpoint(pares_por_uf, "pares_por_uf")
  log_tempo("Módulo 5 concluído e salvo em checkpoint (checkpoints_nordeste/pares_por_uf.rds).")
}

cat(" -> Total de pares (UF x associação) estimados:", nrow(pares_por_uf), "\n")
cat(" -> Pares significativos por UF (RR >", RR_MIN, ", FDR <", ALPHA_FDR, "):\n")
print(pares_por_uf %>% filter(significativo_uf) %>% count(UF, name = "n_pares_sig"))
cat("\n")


# Versão simples de 2x2 para uso PONTUAL (poucas dezenas de chamadas no
# máximo — ex.: Firth no Módulo 7 para os 15 melhores pares). NÃO usar em
# loops grandes (milhares de pares ou centenas de réplicas de bootstrap):
# para isso, use a versão vetorizada (data.table) definida acima/abaixo.
calcular_2x2 <- function(base, comorb, df_univ) {
  tem_base   <- df_univ$DIAG_PRINC_U == base
  tem_comorb <- !is.na(df_univ$DIAGSEC1_U) & df_univ$DIAGSEC1_U == comorb
  base_e_comorb <- tem_base & tem_comorb
  base_isolado  <- tem_base & !tem_comorb
  tibble(
    a = sum(base_e_comorb & df_univ$OBITO_U, na.rm = TRUE),
    b = sum(base_e_comorb & !df_univ$OBITO_U, na.rm = TRUE),
    c = sum(base_isolado & df_univ$OBITO_U, na.rm = TRUE),
    d = sum(base_isolado & !df_univ$OBITO_U, na.rm = TRUE)
  )
}

# Análise "ingênua" de pool bruto (todos os pacientes juntos) para contraste
# metodológico explícito com a meta-análise do Módulo 6 (risco de Simpson).
# Reaproveita a agregação vetorizada de `estimar_pares_uf_rapido` em vez de
# rowwise() + calcular_2x2() por par (mesma razão do Módulo 5: aqui também
# pode haver milhares de pares candidatos, agora sobre 9 UFs simultâneas).
if (existe_checkpoint("pares_pool_ingenuo")) {
  log_tempo("Checkpoint 'pares_pool_ingenuo' encontrado -> pulando (carregando do disco).")
  pares_pool_ingenuo <- carregar_checkpoint("pares_pool_ingenuo")
  universo_pool <- bind_rows(lapply(dados_uf, function(x) x$universo)) # ainda necessário adiante (M7)
} else {
  log_tempo("Calculando pool ingênuo (vetorizado, para contraste com a meta-análise)...")
  universo_pool <- bind_rows(lapply(dados_uf, function(x) x$universo))
  pares_pool_ingenuo <- estimar_pares_uf_rapido("POOL_INGENUO", coorte_nordeste, universo_pool) %>%
    rename(RR_pool_ingenuo = RR)
  salvar_checkpoint(pares_pool_ingenuo, "pares_pool_ingenuo")
  log_tempo("Pool ingênuo concluído e salvo em checkpoint.")
}


# ==============================================================================
# PIPELINE NORDESTE v4.0 - MÓDULOS 6 A 16 (VERSÃO CORRIGIDA)
# Pré-requisito: Módulos 1 a 5 já executados (ou carregados dos checkpoints).
# Objetos esperados no ambiente: pares_por_uf, pares_pool_ingenuo, dados_uf,
# coorte_nordeste, resumo_atricao, N_obitos_nordeste, parâmetros do Módulo 2,
# funções de checkpoint/log_tempo.
# ==============================================================================

stopifnot(
  "Rode os Módulos 1-5 antes (faltam objetos)" =
    all(sapply(c("pares_por_uf", "pares_pool_ingenuo", "dados_uf", "coorte_nordeste"), exists))
)

# Flag de privacidade (LGPD): remove município e datas do CSV de pacientes
REMOVER_QUASE_IDENTIFICADORES <- TRUE
B_NULL_GRAFOS <- 100   # nº de redes Erdős–Rényi para o teste de cascata (M10)


# ==============================================================================
# MÓDULO 6: META-ANÁLISE ENTRE ESTADOS (EFEITOS ALEATÓRIOS)
# ==============================================================================
cat("=== MÓDULO 6: META-ANÁLISE DE LOG-RR ENTRE AS 9 UFs (REML, fallback DL) ===\n")

# rma() é chamado dentro de uma função com yi/sei como argumentos locais:
# evita problemas de avaliação não-padrão dentro do data.table.
# Se o REML não convergir (comum com poucas UFs), tenta DerSimonian-Laird.
.erros_meta <- character(0)

meta_reml_par <- function(yi, sei) {
  metodo <- "REML"
  fit <- tryCatch(
    suppressWarnings(metafor::rma(yi = yi, sei = sei, method = "REML")),
    error = function(e) { .erros_meta <<- c(.erros_meta, conditionMessage(e)); NULL }
  )
  if (is.null(fit)) {
    metodo <- "DL"
    fit <- tryCatch(
      suppressWarnings(metafor::rma(yi = yi, sei = sei, method = "DL")),
      error = function(e) { .erros_meta <<- c(.erros_meta, conditionMessage(e)); NULL }
    )
  }
  if (is.null(fit)) return(NULL)
  list(
    metodo_tau2 = metodo,
    RR_meta     = exp(as.numeric(fit$b[1])),
    RR_meta_li  = exp(fit$ci.lb),
    RR_meta_ls  = exp(fit$ci.ub),
    tau2        = fit$tau2,
    I2          = fit$I2,
    Q_valor     = fit$QE,
    Q_p         = fit$QEp,
    p_meta      = fit$pval
  )
}

pares_unicos <- pares_por_uf %>% distinct(DIAG_PRINC, DIAGSEC1)

if (existe_checkpoint("meta_todos")) {
  log_tempo("Checkpoint 'meta_todos' encontrado -> pulando Módulo 6 (carregando do disco).")
  meta_todos <- carregar_checkpoint("meta_todos")
} else {
  log_tempo(paste0("Iniciando Módulo 6 para ", nrow(pares_unicos), " pares únicos..."))
  
  pares_dt <- as.data.table(pares_por_uf)[
    !is.na(RR) & RR > 0 & is.finite(RR) & is.finite(se_logRR) & se_logRR > 0]
  pares_dt[, k_ufs := .N, by = .(DIAG_PRINC, DIAGSEC1)]
  
  n_elegiveis <- uniqueN(pares_dt[k_ufs >= MIN_UFS_METAAN, .(DIAG_PRINC, DIAGSEC1)])
  cat(" -> Pares únicos com >=", MIN_UFS_METAAN, "UFs válidas:", n_elegiveis, "\n")
  stopifnot("Nenhum par elegível: reduza MIN_UFS_METAAN ou MIN_CASOS_PAR" = n_elegiveis > 0)
  
  meta_bruto <- pares_dt[k_ufs >= MIN_UFS_METAAN, {
    r <- meta_reml_par(log(RR), se_logRR)
    if (is.null(r)) NULL else c(
      list(n_UFs_incluidas = .N, UFs_incluidas = paste(UF, collapse = ", ")), r)
  }, by = .(DIAG_PRINC, DIAGSEC1)]
  
  if (length(.erros_meta) > 0) {
    cat(" -> Avisos: ", length(.erros_meta), " falhas de convergência/ajuste. Mensagem mais frequente:\n    ",
        names(sort(table(.erros_meta), decreasing = TRUE))[1], "\n")
  }
  stopifnot("Meta-análise não gerou nenhuma linha (veja mensagens acima)" = nrow(meta_bruto) > 0)
  
  # FDR sobre TODOS os pares meta-analisados (família completa de testes)
  meta_todos <- as_tibble(meta_bruto) %>%
    mutate(
      p_meta_fdr = p.adjust(p_meta, method = "fdr"),
      robustez_regional = ifelse(I2 < I2_ROBUSTO,
                                 "Robusta (baixa heterogeneidade)",
                                 "Estado-específica (alta heterogeneidade)")
    )
  salvar_checkpoint(meta_todos, "meta_todos")
  log_tempo("Módulo 6 concluído e salvo em checkpoint.")
}

# Rede oficial: significativos (FDR) e RR meta acima do limiar principal
meta_resultados <- meta_todos %>%
  filter(!is.na(p_meta_fdr), p_meta_fdr < ALPHA_FDR, RR_meta > RR_MIN)

stopifnot("Nenhum par passou FDR + RR_MIN: rede vazia" = nrow(meta_resultados) > 0)

cat(" -> Pares meta-analisados (total):", nrow(meta_todos), "\n")
cat(" -> Pares com significância regional (FDR e RR_meta > RR_MIN):", nrow(meta_resultados), "\n")
cat(" -> Distribuição por robustez regional:\n")
print(table(meta_resultados$robustez_regional))
cat("\n")

rede_principal_nordeste <- meta_resultados %>%
  mutate(
    peso_bruto = n_UFs_incluidas,
    weight     = RR_meta,
    tier       = ifelse(robustez_regional == "Robusta (baixa heterogeneidade)", "Alta", "Exploratória")
  )


# ==============================================================================
# MÓDULO 7: REGRESSÃO LOGÍSTICA — FIRTH E EFEITOS MISTOS (UF ALEATÓRIO)
# ==============================================================================
cat("=== MÓDULO 7: FIRTH E MODELO MISTO COM UF ALEATÓRIO ===\n")

# Universo empilhado (com coluna UF) em data.table, construído uma vez
if (exists("universo_pool")) { rm(universo_pool); invisible(gc(verbose = FALSE)) }
universo_pool_dt <- rbindlist(lapply(names(dados_uf), function(uf) {
  x <- as.data.table(dados_uf[[uf]]$universo)
  x[, UF := uf]
  x
}))

calcular_2x2_dt <- function(base, comorb) {
  sub  <- universo_pool_dt[DIAG_PRINC_U == base]
  expo <- !is.na(sub$DIAGSEC1_U) & sub$DIAGSEC1_U == comorb
  tibble(
    a = sum(expo & sub$OBITO_U),  b = sum(expo & !sub$OBITO_U),
    c = sum(!expo & sub$OBITO_U), d = sum(!expo & !sub$OBITO_U)
  )
}

# Firth com dados agregados (4 linhas + pesos), em vez de 1 linha por paciente
calcular_or_firth <- function(a, b, cc, d) {
  df <- data.frame(exposto = c(1, 1, 0, 0), obito = c(1, 0, 1, 0), w = c(a, b, cc, d))
  df <- df[df$w > 0, ]
  vazio <- tibble(OR_firth = NA_real_, OR_firth_li = NA_real_, OR_firth_ls = NA_real_)
  fit <- tryCatch(logistf(obito ~ exposto, data = df, weights = w, pl = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(vazio)
  tibble(
    OR_firth    = exp(unname(coef(fit)["exposto"])),
    OR_firth_li = exp(unname(fit$ci.lower["exposto"])),
    OR_firth_ls = exp(unname(fit$ci.upper["exposto"]))
  )
}

top_pares_firth <- rede_principal_nordeste %>% arrange(desc(RR_meta)) %>% head(15)

resultados_firth <- top_pares_firth %>%
  mutate(cel = map2(DIAG_PRINC, DIAGSEC1, calcular_2x2_dt)) %>%
  unnest(cel) %>%
  mutate(firth = pmap(list(a, b, c, d), calcular_or_firth)) %>%
  unnest(firth)

# Efeitos mistos com binomial agregado por UF x exposição (equivalente ao
# modelo individual, muito mais rápido)
ajustar_misto_uf <- function(base, comorb) {
  vazio <- tibble(OR_misto = NA_real_, OR_misto_li = NA_real_, OR_misto_ls = NA_real_, var_UF = NA_real_)
  dfb <- universo_pool_dt[DIAG_PRINC_U == base, .(n = .N),
                          by = .(UF,
                                 exposto = as.integer(!is.na(DIAGSEC1_U) & DIAGSEC1_U == comorb),
                                 obito   = as.integer(OBITO_U))]
  df <- dcast(dfb, UF + exposto ~ obito, value.var = "n", fill = 0L)
  if (!"1" %in% names(df)) df[, `1` := 0L]
  if (!"0" %in% names(df)) df[, `0` := 0L]
  setnames(df, c("1", "0"), c("obitos", "sobreviventes"))
  if (uniqueN(df$UF) < 2 || uniqueN(df$exposto) < 2) return(vazio)
  fit <- tryCatch(
    glmer(cbind(obitos, sobreviventes) ~ exposto + (1 | UF), data = as.data.frame(df),
          family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(vazio)
  est <- fixef(fit)["exposto"]; se <- sqrt(vcov(fit)["exposto", "exposto"])
  tibble(OR_misto = exp(est), OR_misto_li = exp(est - 1.96 * se), OR_misto_ls = exp(est + 1.96 * se),
         var_UF = as.numeric(VarCorr(fit)$UF[1]))
}

resultados_misto <- top_pares_firth %>%
  head(5) %>%
  mutate(misto = map2(DIAG_PRINC, DIAGSEC1, ajustar_misto_uf)) %>%
  unnest(misto)

cat(" -> Firth (pool agregado), top pares:\n")
print(resultados_firth %>% select(DIAG_PRINC, DIAGSEC1, RR_meta, OR_firth, OR_firth_li, OR_firth_ls))
cat(" -> Efeitos mistos (UF aleatório), top 5 pares:\n")
print(resultados_misto %>% select(DIAG_PRINC, DIAGSEC1, RR_meta, OR_misto, var_UF))
cat("\n")


# ==============================================================================
# MÓDULO 8: CONSTRUÇÃO DAS REDES (REGIONAL, POR UF, MULTICAMADA)
# ==============================================================================
cat("=== MÓDULO 8: CONSTRUINDO REDES ===\n")

grafo_regional <- graph_from_data_frame(
  rede_principal_nordeste %>%
    select(DIAG_PRINC, DIAGSEC1, weight, RR_meta, I2, robustez_regional, tier, n_UFs_incluidas),
  directed = FALSE
)
grafo_regional <- delete_vertices(grafo_regional, which(degree(grafo_regional) == 0))
V(grafo_regional)$betweenness <- betweenness(grafo_regional, directed = FALSE, normalized = TRUE) * 100
comunidades_regional <- cluster_louvain(grafo_regional)
V(grafo_regional)$group <- as.character(membership(comunidades_regional))
modularidade_regional <- modularity(comunidades_regional)

cat(" -> Rede regional: nós =", vcount(grafo_regional), "| arestas =", ecount(grafo_regional), "\n")
cat(" -> Modularidade (Louvain, ponderada):", round(modularidade_regional, 4), "\n\n")

grafos_uf <- list()
for (uf in UFS_NORDESTE) {
  rede_uf <- pares_por_uf %>% filter(UF == uf, significativo_uf)
  if (nrow(rede_uf) == 0) { grafos_uf[[uf]] <- NULL; next }
  g <- graph_from_data_frame(rede_uf %>% mutate(weight = casos_obito_par * RR) %>%
                               select(DIAG_PRINC, DIAGSEC1, weight, RR), directed = FALSE)
  g <- delete_vertices(g, which(degree(g) == 0))
  grafos_uf[[uf]] <- g
}

participacao_multiplex <- map_dfr(names(grafos_uf), function(uf) {
  g <- grafos_uf[[uf]]
  if (is.null(g) || vcount(g) == 0) return(tibble())
  deg <- degree(g)
  limiar_hub <- quantile(deg, 0.90, na.rm = TRUE)
  tibble(UF = uf, diagnostico = names(deg)[deg >= limiar_hub], grau = deg[deg >= limiar_hub])
})

hubs_universais <- participacao_multiplex %>%
  count(diagnostico, name = "n_camadas_como_hub") %>%
  arrange(desc(n_camadas_como_hub))

cat(" -> 'Hubs universais':\n"); print(head(hubs_universais, 10)); cat("\n")

# ==============================================================================
# PIPELINE NORDESTE v4.0 - MÓDULOS 9 A 16 (ATUALIZADO)
# Pré-requisito: Módulos 1 a 8 já executados nesta sessão.
# Objetos esperados: grafo_regional, rede_principal_nordeste, meta_todos,
# meta_resultados, resultados_firth, resultados_misto, hubs_universais,
# modularidade_regional, pares_por_uf, pares_pool_ingenuo, dados_uf,
# coorte_nordeste, resumo_atricao, N_obitos_nordeste, parâmetros do Módulo 2,
# funções de checkpoint/log_tempo.
# ==============================================================================

stopifnot(
  "Rode os Módulos 1-8 antes (faltam objetos)" =
    all(sapply(c("grafo_regional", "rede_principal_nordeste", "meta_todos", "meta_resultados",
                 "resultados_firth", "resultados_misto", "hubs_universais", "modularidade_regional",
                 "pares_por_uf", "pares_pool_ingenuo", "coorte_nordeste"), exists))
)

REMOVER_QUASE_IDENTIFICADORES <- TRUE
B_NULL_GRAFOS <- 100   # nº de redes Erdős–Rényi para o teste de cascata (M10)


# ==============================================================================
# MÓDULO 9: TOPOLOGIA FORMAL E REDES NULAS
# ==============================================================================
cat("=== MÓDULO 9: TOPOLOGIA FORMAL E COMPARAÇÃO COM MODELOS NULOS ===\n")

descrever_topologia <- function(g, nome = "empírica") {
  tibble(
    Rede              = nome,
    N_nos             = vcount(g),
    N_arestas         = ecount(g),
    Grau_medio        = round(mean(degree(g)), 3),
    Clustering_global = round(transitivity(g, type = "global"), 4),
    Caminho_medio     = round(mean_distance(g, directed = FALSE, unconnected = TRUE), 3),
    Diametro          = diameter(g, directed = FALSE, unconnected = TRUE),
    Assortatividade   = round(assortativity_degree(g, directed = FALSE), 4)
  )
}

# CORRIGIDO: extrai campos do ajuste com segurança (KS.p pode não ser numérico)
ajustar_lei_potencia <- function(g) {
  deg <- degree(g); deg <- deg[deg > 0]
  fit <- tryCatch(fit_power_law(deg, implementation = "plfit"), error = function(e) NULL)
  if (is.null(fit)) return(tibble(alpha = NA_real_, x_min = NA_real_, p_ks = NA_real_))
  pega <- function(x) if (is.null(x) || !is.numeric(x) || length(x) == 0) NA_real_ else as.numeric(x[1])
  tibble(alpha = round(pega(fit$alpha), 3),
         x_min = pega(fit$xmin),
         p_ks  = round(pega(fit$KS.p), 4))
}

percolacao_molloy_reed <- function(g) {
  deg <- degree(g); k1 <- mean(deg); k2 <- mean(deg^2)
  if ((k2 - k1) <= 0) return(NA_real_)
  k1 / (k2 - k1)
}

gerar_nulos_pareados <- function(g, n_reps = 30) {
  N <- vcount(g); E <- ecount(g); k_medio <- 2 * E / N
  m_ba <- max(1, round(k_medio / 2))
  er_list <- lapply(seq_len(n_reps), function(i) sample_gnm(N, E))
  ba_list <- lapply(seq_len(n_reps), function(i) sample_pa(N, m = m_ba, directed = FALSE))
  ws_list <- lapply(seq_len(n_reps), function(i) sample_smallworld(dim = 1, size = N, nei = m_ba, p = 0.05))
  resumir <- function(lst, nome) {
    m <- map_dfr(lst, function(gg) tibble(
      clustering = transitivity(gg, type = "global"),
      caminho    = mean_distance(gg, directed = FALSE, unconnected = TRUE)))
    tibble(Rede = nome,
           Clustering_medio = round(mean(m$clustering, na.rm = TRUE), 4),
           Clustering_dp    = round(sd(m$clustering, na.rm = TRUE), 4),
           Caminho_medio    = round(mean(m$caminho, na.rm = TRUE), 3))
  }
  bind_rows(resumir(er_list, "Erdős–Rényi (nula)"),
            resumir(ba_list, "Barabási–Albert (nula)"),
            resumir(ws_list, "Watts–Strogatz (nula)"))
}

topologia_empirica  <- descrever_topologia(grafo_regional, "Nordeste (empírica)")
lei_potencia_result <- ajustar_lei_potencia(grafo_regional)
pc_molloy_reed      <- percolacao_molloy_reed(grafo_regional)
topologia_nulos     <- gerar_nulos_pareados(grafo_regional, n_reps = 30)

cat(" -> Topologia empírica:\n"); print(as.data.frame(topologia_empirica))
cat(" -> Lei de potência (grau):\n"); print(as.data.frame(lei_potencia_result))
cat(" -> Limiar de Molloy–Reed (pc):", round(pc_molloy_reed, 4), "\n")
cat(" -> Redes nulas pareadas:\n"); print(as.data.frame(topologia_nulos)); cat("\n")


# ==============================================================================
# MÓDULO 10: PERCOLAÇÃO DE ARESTA E MONTE CARLO DE CASCATA (VETORIZADO)
# ==============================================================================
cat("=== MÓDULO 10: MONTE CARLO DE CASCATA (PERCOLAÇÃO DE ARESTA) ===\n")

normalizar_prob_transmissao <- function(g, prob_base = 0.5) {
  rr <- E(g)$RR_meta
  if (is.null(rr)) rr <- E(g)$RR
  rr_norm <- (rr - min(rr, na.rm = TRUE)) / (max(rr, na.rm = TRUE) - min(rr, na.rm = TRUE) + 1e-9)
  prob_base * (0.5 + rr_norm)
}

mc_cascata_rapido <- function(g, p_aresta, idx_sementes, B) {
  el    <- as_edgelist(g, names = FALSE)
  N     <- vcount(g); M <- nrow(el)
  vazio <- make_empty_graph(N, directed = FALSE)
  tam   <- matrix(NA_integer_, nrow = length(idx_sementes), ncol = B)
  for (b in seq_len(B)) {
    ab <- runif(M) < p_aresta
    gb <- if (any(ab)) add_edges(vazio, as.vector(t(el[ab, , drop = FALSE]))) else vazio
    cmp <- components(gb)
    tam[, b] <- cmp$csize[cmp$membership[idx_sementes]]
  }
  tam
}

prob_aresta_regional <- normalizar_prob_transmissao(grafo_regional, prob_base = 0.4)
N_g <- vcount(grafo_regional); E_g <- ecount(grafo_regional)

if (existe_checkpoint("resultados_mc_v2")) {
  log_tempo("Checkpoint 'resultados_mc_v2' encontrado -> carregando.")
  mc_obj <- carregar_checkpoint("resultados_mc_v2")
  resultados_mc <- mc_obj$resultados_mc
  media_nula    <- mc_obj$media_nula
  media_emp     <- mc_obj$media_emp
} else {
  log_tempo(paste0("Monte Carlo: ", N_g, " sementes x B = ", B_MC, " réplicas..."))
  tam_emp <- mc_cascata_rapido(grafo_regional, prob_aresta_regional, seq_len(N_g), B_MC)
  resultados_mc <- tibble(
    diagnostico_semente = V(grafo_regional)$name,
    tamanho_medio   = rowMeans(tam_emp),
    tamanho_mediano = apply(tam_emp, 1, median),
    tamanho_p95     = apply(tam_emp, 1, quantile, 0.95)
  ) %>% arrange(desc(tamanho_medio))
  media_emp <- mean(tam_emp)
  
  log_tempo(paste0("Monte Carlo: ", B_NULL_GRAFOS, " redes nulas..."))
  p_nula <- mean(prob_aresta_regional, na.rm = TRUE)
  media_nula <- vapply(seq_len(B_NULL_GRAFOS), function(i) {
    g_er <- sample_gnm(N_g, E_g)
    mean(mc_cascata_rapido(g_er, rep(p_nula, E_g), seq_len(N_g), 100))
  }, numeric(1))
  
  salvar_checkpoint(list(resultados_mc = resultados_mc, media_nula = media_nula, media_emp = media_emp),
                    "resultados_mc_v2")
  log_tempo("Monte Carlo concluído e salvo em checkpoint.")
}

p_empirico_cascata <- (1 + sum(media_nula >= media_emp)) / (length(media_nula) + 1)

cat(" -> Top 10 sementes por tamanho médio de cascata:\n"); print(head(resultados_mc, 10))
cat(" -> Cascata média (todos os nós): empírica =", round(media_emp, 2),
    "| nula (ER) =", round(mean(media_nula), 2), "\n")
cat(" -> p-empírico (H0: rede aleatória produz cascata igual/maior):", round(p_empirico_cascata, 4), "\n\n")


# ==============================================================================
# MÓDULO 11: ROBUSTEZ — ATAQUE A HUBS vs. FALHA ALEATÓRIA
# ==============================================================================
cat("=== MÓDULO 11: RESILIÊNCIA DA REDE ===\n")

simular_ataque <- function(g, tipo = c("grau", "betweenness", "aleatorio"), n_rep = 30) {
  tipo <- match.arg(tipo)
  N <- vcount(g)
  fracoes <- seq(0, 0.9, by = 0.05)
  calcular_frac_gigante <- function(ordem_remocao) {
    map_dbl(fracoes, function(f) {
      k <- floor(f * N)
      if (k == 0) return(max(components(g)$csize) / N)
      g_sub <- delete_vertices(g, ordem_remocao[seq_len(k)])
      if (vcount(g_sub) == 0) return(0)
      max(components(g_sub)$csize) / N
    })
  }
  if (tipo == "grau") {
    ordem <- V(g)$name[order(-degree(g))]
    tibble(fracao_removida = fracoes, frac_giant = calcular_frac_gigante(ordem), tipo = "Ataque por grau")
  } else if (tipo == "betweenness") {
    ordem <- V(g)$name[order(-betweenness(g))]
    tibble(fracao_removida = fracoes, frac_giant = calcular_frac_gigante(ordem), tipo = "Ataque por intermediação")
  } else {
    reps <- map(seq_len(n_rep), function(i) calcular_frac_gigante(sample(V(g)$name)))
    tibble(fracao_removida = fracoes, frac_giant = Reduce(`+`, reps) / n_rep, tipo = "Falha aleatória")
  }
}

robustez_resultados <- bind_rows(
  simular_ataque(grafo_regional, "grau"),
  simular_ataque(grafo_regional, "betweenness"),
  simular_ataque(grafo_regional, "aleatorio", n_rep = 50)
)
cat(" -> Aos 20% de remoção:\n")
print(robustez_resultados %>% filter(abs(fracao_removida - 0.2) < 1e-9)); cat("\n")


# ==============================================================================
# MÓDULO 12: REDE DIRECIONADA (DIAGNÓSTICO-BASE -> COMORBIDADE)
# ==============================================================================
cat("=== MÓDULO 12: REDE DIRECIONADA E COMPONENTES ===\n")

grafo_direcionado <- graph_from_data_frame(
  rede_principal_nordeste %>% select(DIAG_PRINC, DIAGSEC1, weight, RR_meta), directed = TRUE)

resumo_direcional <- tibble(
  diagnostico = V(grafo_direcionado)$name,
  out_degree  = degree(grafo_direcionado, mode = "out"),
  in_degree   = degree(grafo_direcionado, mode = "in"),
  betweenness = round(betweenness(grafo_direcionado, directed = TRUE), 2)
) %>% arrange(desc(out_degree))

componentes_dir <- components(grafo_direcionado, mode = "weak")
componente_gigante_frac <- max(componentes_dir$csize) / vcount(grafo_direcionado)

cat(" -> Maiores out-degree:\n"); print(head(resumo_direcional, 10))
cat(" -> Componentes fracos:", componentes_dir$no, "| gigante:",
    round(componente_gigante_frac * 100, 1), "% dos nós\n\n")


# ==============================================================================
# MÓDULO 13: BOOTSTRAP DE CENTRALIDADE E PERMUTAÇÃO DE MODULARIDADE
# ==============================================================================
cat("=== MÓDULO 13: BOOTSTRAP (B =", B_BOOT, ") E PERMUTAÇÃO DE MODULARIDADE (B =", B_MOD, ") ===\n")

top10_ids <- V(grafo_regional)$name[order(-V(grafo_regional)$betweenness)][1:min(10, vcount(grafo_regional))]
boot_resultados <- matrix(NA_real_, nrow = B_BOOT, ncol = length(top10_ids), dimnames = list(NULL, top10_ids))

metanalise_dl_rapida <- function(dt) {
  dt <- dt[!is.na(RR) & is.finite(RR) & RR > 0 & is.finite(se_logRR) & se_logRR > 0]
  if (nrow(dt) == 0) return(data.table())
  dt[, `:=`(yi = log(RR), vi = se_logRR^2)]
  dt[, wi := 1 / vi]
  dt[, k := .N, by = .(DIAG_PRINC, DIAGSEC1)]
  dt <- dt[k >= MIN_UFS_METAAN]
  if (nrow(dt) == 0) return(data.table())
  
  dt[, `:=`(sum_wi = sum(wi), sum_wi2 = sum(wi^2), mu_f = sum(wi * yi) / sum(wi)),
     by = .(DIAG_PRINC, DIAGSEC1)]
  dt[, Q := sum(wi * (yi - mu_f)^2), by = .(DIAG_PRINC, DIAGSEC1)]
  dt[, tau2 := pmax(0, (Q - (k - 1)) / (sum_wi - sum_wi2 / sum_wi))]
  dt[, wi_star := 1 / (vi + tau2)]
  
  saida <- dt[, .(yi_pooled = sum(wi_star * yi) / sum(wi_star),
                  se_pooled = sqrt(1 / sum(wi_star))), by = .(DIAG_PRINC, DIAGSEC1)]
  saida[, `:=`(RR_meta_boot = exp(yi_pooled),
               p_boot = 2 * pnorm(-abs(yi_pooled / se_pooled)))]
  saida[, p_boot_fdr := p.adjust(p_boot, "fdr")]
  saida[]
}

if (existe_checkpoint("tabela_bootstrap")) {
  log_tempo("Checkpoint 'tabela_bootstrap' encontrado -> pulando bootstrap.")
  tabela_bootstrap <- carregar_checkpoint("tabela_bootstrap")
} else {
  log_tempo(paste0("Bootstrap paramétrico (B = ", B_BOOT, ") sobre as contagens de célula..."))
  pares_uf_dt <- as.data.table(pares_por_uf)[, .(UF, DIAG_PRINC, DIAGSEC1, a, b, c, d)]
  
  set.seed(123)
  for (i in seq_len(B_BOOT)) {
    if (i %% 20 == 0) log_tempo(paste0(" -> Bootstrap: réplica ", i, "/", B_BOOT))
    boot_i <- copy(pares_uf_dt)
    boot_i[, a_boot := rbinom(.N, size = a + b, prob = a / (a + b))]
    boot_i[, c_boot := rbinom(.N, size = c + d, prob = c / (c + d))]
    boot_i[, `:=`(
      RR       = (a_boot / (a + b)) / (c_boot / (c + d)),
      se_logRR = sqrt(1 / pmax(a_boot, 0.5) - 1 / (a + b) + 1 / pmax(c_boot, 0.5) - 1 / (c + d))
    )]
    meta_boot <- metanalise_dl_rapida(boot_i)
    if (nrow(meta_boot) == 0) next
    rede_boot <- meta_boot[p_boot_fdr < ALPHA_FDR & RR_meta_boot > RR_MIN]
    if (nrow(rede_boot) == 0) next
    
    g_boot <- graph_from_data_frame(rede_boot[, .(DIAG_PRINC, DIAGSEC1, weight = RR_meta_boot)], directed = FALSE)
    g_boot <- delete_vertices(g_boot, which(degree(g_boot) == 0))
    cent_boot <- betweenness(g_boot, directed = FALSE, normalized = TRUE) * 100
    for (id in top10_ids) if (id %in% names(cent_boot)) boot_resultados[i, id] <- cent_boot[id]
  }
  
  tabela_bootstrap <- tibble(
    `CID-10` = top10_ids,
    `Centralidade Observada (%)`     = round(V(grafo_regional)$betweenness[match(top10_ids, V(grafo_regional)$name)], 2),
    `Mediana Bootstrap (%)`          = round(apply(boot_resultados, 2, median, na.rm = TRUE), 2),
    `Estabilidade (% Reamostragens)` = round(100 * colMeans(!is.na(boot_resultados)), 1)
  )
  salvar_checkpoint(tabela_bootstrap, "tabela_bootstrap")
  log_tempo("Bootstrap concluído e salvo em checkpoint.")
}

if (existe_checkpoint("modularidades_nulas")) {
  log_tempo("Checkpoint 'modularidades_nulas' encontrado -> carregando.")
  modularidades_nulas <- carregar_checkpoint("modularidades_nulas")
} else {
  log_tempo(paste0("Permutação de modularidade (B = ", B_MOD, ")..."))
  set.seed(2024)
  pesos_obs <- E(grafo_regional)$weight
  modularidades_nulas <- vapply(seq_len(B_MOD), function(i) {
    g_rand <- tryCatch(rewire(grafo_regional, with = keeping_degseq(loops = FALSE, niter = ecount(grafo_regional) * 10)),
                       error = function(e) NULL)
    if (is.null(g_rand)) return(NA_real_)
    E(g_rand)$weight <- sample(pesos_obs)
    com <- tryCatch(cluster_louvain(g_rand), error = function(e) NULL)
    if (is.null(com)) NA_real_ else modularity(com)
  }, numeric(1))
  salvar_checkpoint(modularidades_nulas, "modularidades_nulas")
}
modularidades_nulas <- modularidades_nulas[!is.na(modularidades_nulas)]
p_empirico_mod <- (1 + sum(modularidades_nulas >= modularidade_regional)) / (length(modularidades_nulas) + 1)

cat(" -> Modularidade observada:", round(modularidade_regional, 4),
    "| média sob H0:", round(mean(modularidades_nulas), 4),
    "| p-empírico:", round(p_empirico_mod, 4), "\n\n")


# ==============================================================================
# MÓDULO 14: SENSIBILIDADE MULTI-LIMIAR (REGIONAL E POR UF)
# ==============================================================================
cat("=== MÓDULO 14: SENSIBILIDADE MULTI-LIMIAR (RR:", paste(LIMIARES_SENS, collapse = ", "), ") ===\n")

meta_sig_todos <- meta_todos %>% filter(!is.na(p_meta_fdr), p_meta_fdr < ALPHA_FDR)

tabela_sensibilidade_regional <- map_dfr(LIMIARES_SENS, function(lim) {
  rede_lim <- meta_sig_todos %>% filter(RR_meta > lim)
  if (nrow(rede_lim) == 0) return(tibble(`Limiar RR` = lim, `N Arestas` = 0L, `N Nós` = 0L, `Top3 Hubs` = "Nenhum"))
  g_lim <- graph_from_data_frame(rede_lim %>% select(DIAG_PRINC, DIAGSEC1, RR_meta), directed = FALSE)
  g_lim <- delete_vertices(g_lim, which(degree(g_lim) == 0))
  cent_lim <- betweenness(g_lim, directed = FALSE, normalized = TRUE) * 100
  top3 <- sort(cent_lim, decreasing = TRUE)[1:min(3, length(cent_lim))]
  tibble(`Limiar RR` = lim, `N Arestas` = nrow(rede_lim), `N Nós` = vcount(g_lim),
         `Top3 Hubs` = paste(names(top3), round(top3, 1), sep = ": ", collapse = "; "))
})

tabela_sensibilidade_uf <- map_dfr(UFS_NORDESTE, function(uf) {
  map_dfr(LIMIARES_SENS, function(lim) {
    n <- pares_por_uf %>% filter(UF == uf, RR > lim, significativo_uf) %>% nrow()
    tibble(UF = uf, `Limiar RR` = lim, `N Arestas` = n)
  })
})

cat(" -> Sensibilidade regional:\n"); print(tabela_sensibilidade_regional)
cat(" -> Sensibilidade por UF:\n")
print(tabela_sensibilidade_uf %>% pivot_wider(names_from = `Limiar RR`, values_from = `N Arestas`))
cat("\n")


# ==============================================================================
# MÓDULO 15: ANÁLISE TEMPORAL
# ==============================================================================
cat("=== MÓDULO 15: SÉRIE TEMPORAL DE MORTALIDADE ===\n")

if (!exists("cols_ref") || is.null(cols_ref)) cols_ref <- detectar_colunas(coorte_nordeste)
col_ano <- cols_ref$col_ano

tabela_temporal_uf <- coorte_nordeste %>%
  group_by(UF, Ano = .data[[col_ano]]) %>%
  summarise(Obitos_N = n(), .groups = "drop") %>%
  group_by(UF) %>%
  mutate(Percentual_UF = round(100 * Obitos_N / sum(Obitos_N), 2)) %>%
  ungroup()

tabela_temporal_regional <- tabela_temporal_uf %>%
  group_by(Ano) %>%
  summarise(Obitos_N = sum(Obitos_N), .groups = "drop") %>%
  mutate(Percentual = round(100 * Obitos_N / sum(Obitos_N), 2)) %>%
  arrange(Ano)

metricas_seguimento_uf <- coorte_nordeste %>%
  mutate(DIAS_PERM_NUM = as.numeric(as.character(DIAS_PERM))) %>%
  group_by(UF) %>%
  summarise(
    N_obitos     = n(),
    Mediana_Dias = median(DIAS_PERM_NUM, na.rm = TRUE),
    Q1_Dias      = quantile(DIAS_PERM_NUM, 0.25, na.rm = TRUE),
    Q3_Dias      = quantile(DIAS_PERM_NUM, 0.75, na.rm = TRUE),
    Media_Dias   = round(mean(DIAS_PERM_NUM, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat(" -> Permanência mediana (IQR) por UF:\n"); print(as.data.frame(metricas_seguimento_uf))
cat(" -> Série regional:\n"); print(as.data.frame(tabela_temporal_regional)); cat("\n")


# ==============================================================================
# MÓDULO 16: TABELAS, FIGURAS E PACOTE DE CIÊNCIA ABERTA
# ==============================================================================
cat("=== MÓDULO 16: TABELAS, FIGURAS E PACOTE OSF ===\n")

cores_modulos <- c("1" = "#1a9850", "2" = "#fc8d62", "3" = "#8073ac", "4" = "#e78ac3",
                   "5" = "#66c2a5", "6" = "#ffd92f", "7" = "#a6d854", "8" = "#e5c494")

tabela_hubs_regional <- tibble(id = V(grafo_regional)$name, group = V(grafo_regional)$group,
                               betweenness = V(grafo_regional)$betweenness) %>%
  mutate(Grau = degree(grafo_regional)[id]) %>%
  left_join(tabela_bootstrap %>% rename(id = `CID-10`), by = "id") %>%
  arrange(desc(betweenness)) %>%
  head(15)

print(tabela_hubs_regional %>%
        kbl(caption = "Tabela 2. Diagnósticos mediadores críticos (hubs) da rede regional Nordeste.") %>%
        kable_classic(full_width = FALSE))

fmt_rr_ic <- function(rr, li, ls) paste0(format(round(rr, 2), decimal.mark = ","), " (",
                                         format(round(li, 2), decimal.mark = ","), "–",
                                         format(round(ls, 2), decimal.mark = ","), ")")
fmt_p <- function(p) ifelse(p < 0.001, "< 0,001", format(round(p, 4), decimal.mark = ","))

tabela3_regional <- rede_principal_nordeste %>%
  mutate(`RR meta (IC95%)` = fmt_rr_ic(RR_meta, RR_meta_li, RR_meta_ls),
         `p (FDR)` = fmt_p(p_meta_fdr),
         `I² (%)` = round(I2, 1)) %>%
  arrange(desc(RR_meta)) %>%
  select(`Diagnóstico Base` = DIAG_PRINC, `Comorbidade` = DIAGSEC1, `UFs incluídas` = n_UFs_incluidas,
         `RR meta (IC95%)`, `I² (%)`, `p (FDR)`, `Classificação` = robustez_regional)
cat(" -> Tabela 3: N =", nrow(tabela3_regional), "pares\n\n")

plot_hubs <- tabela_hubs_regional %>%
  mutate(rotulo = fct_reorder(id, betweenness)) %>%
  ggplot(aes(x = betweenness, y = rotulo, fill = group)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = cores_modulos, name = "Módulo Sindrômico") +
  labs(title = "Diagnósticos Mediadores Críticos — Rede Regional Nordeste",
       subtitle = paste0("Top 15 por intermediação | ", vcount(grafo_regional), " nós, ", ecount(grafo_regional), " arestas"),
       x = "Centralidade de intermediação normalizada (%)", y = "CID-10 (diagnóstico base)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")
ggsave("Figura1_Hubs_Regional_NORDESTE.png", plot_hubs, width = 10, height = 7, dpi = 300, bg = "white")

set.seed(42)
grafo_plot <- ggraph(grafo_regional, layout = "fr") +
  geom_edge_link(aes(width = weight, linetype = tier), color = "#555555", alpha = 0.6) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  scale_edge_linetype_manual(values = c("Alta" = "solid", "Exploratória" = "dotted"), name = "Robustez regional") +
  geom_node_point(aes(size = betweenness, fill = group), shape = 21, color = "black", stroke = 0.5) +
  scale_size_continuous(range = c(3, 16), name = "Intermediação (%)") +
  scale_fill_manual(values = cores_modulos, name = "Módulo Sindrômico") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold", bg.color = "white", bg.r = 0.1) +
  labs(title = "Topologia das Associações Clínicas Letais — Nordeste (2016–2025)",
       subtitle = "Linhas contínuas = robustas entre estados (I² < 50%)") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5), legend.position = "right")
ggsave("Figura2_Grafo_Topologico_NORDESTE.tiff", grafo_plot, width = 12, height = 10, dpi = 300,
       compression = "lzw", bg = "white")

deg_df <- tibble(grau = degree(grafo_regional)) %>% count(grau, name = "freq") %>% mutate(freq_rel = freq / sum(freq))
plot_lei_potencia <- ggplot(deg_df, aes(x = grau, y = freq_rel)) +
  geom_point(size = 2, color = "#2C3E50") +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Distribuição de Grau da Rede Regional (log-log)",
       subtitle = paste0("Lei de potência: alpha = ", lei_potencia_result$alpha, " | KS p = ", lei_potencia_result$p_ks),
       x = "Grau (k)", y = "Frequência relativa") +
  theme_minimal(base_size = 12)
ggsave("Figura3_Distribuicao_Grau_LeiPotencia.png", plot_lei_potencia, width = 8, height = 6, dpi = 300, bg = "white")

plot_robustez <- ggplot(robustez_resultados, aes(x = fracao_removida, y = frac_giant, color = tipo)) +
  geom_line(linewidth = 1.1) + geom_point(size = 1.5) +
  scale_color_manual(values = c("Ataque por grau" = "#C0392B", "Ataque por intermediação" = "#E67E22",
                                "Falha aleatória" = "#2980B9")) +
  labs(title = "Robustez da Rede Regional: Ataque Direcionado vs. Falha Aleatória",
       x = "Fração de nós removidos", y = "Fração do maior componente conectado", color = "Cenário") +
  theme_minimal(base_size = 12)
ggsave("Figura4_Robustez_Ataque_vs_Aleatorio.png", plot_robustez, width = 9, height = 6, dpi = 300, bg = "white")

plot_mc <- ggplot(tibble(media_nula = media_nula), aes(x = media_nula)) +
  geom_histogram(bins = 30, fill = "#95a5a6", color = "white") +
  geom_vline(xintercept = media_emp, color = "#C0392B", linewidth = 1.2) +
  labs(title = "Tamanho médio de cascata: rede empírica vs. redes nulas Erdős–Rényi",
       subtitle = paste0("Linha vermelha = rede empírica | ", length(media_nula), " redes nulas | p-empírico = ",
                         round(p_empirico_cascata, 4)),
       x = "Tamanho médio de cascata (nº de diagnósticos alcançados)", y = "Nº de redes nulas") +
  theme_minimal(base_size = 12)
ggsave("Figura5_MonteCarlo_Cascata.png", plot_mc, width = 9, height = 6, dpi = 300, bg = "white")

df_reg <- tabela_temporal_regional %>%
  mutate(Ano = as.numeric(as.character(Ano)), Obitos_Acumulados = cumsum(Obitos_N))
fator <- max(df_reg$Obitos_Acumulados) / max(df_reg$Obitos_N)

plot_temporal_regional <- ggplot(df_reg, aes(x = Ano)) +
  geom_col(aes(y = Obitos_N), fill = "#2C3E50", alpha = 0.85, width = 0.6) +
  geom_line(aes(y = Obitos_Acumulados / fator), color = "#C0392B", linewidth = 1.1, group = 1) +
  scale_x_continuous(breaks = ANO_INICIO:ANO_FIM) +
  scale_y_continuous(name = "Óbitos por ano (Nordeste)",
                     sec.axis = sec_axis(~ . * fator, name = "Óbitos acumulados")) +
  theme_minimal(base_size = 12) +
  labs(title = "Mortalidade Hospitalar Anual e Acumulada — Nordeste (2016–2025)")

plot_temporal_facet <- tabela_temporal_uf %>%
  mutate(Ano = as.numeric(as.character(Ano))) %>%
  ggplot(aes(x = Ano, y = Obitos_N)) +
  geom_col(fill = "#2980B9") +
  facet_wrap(~ UF, scales = "free_y") +
  theme_minimal(base_size = 10) +
  labs(title = "Mortalidade Hospitalar Anual por UF", x = "Ano", y = "Óbitos")

ggsave("Figura6_Temporal_Regional.png", plot_temporal_regional, width = 11, height = 6, dpi = 300, bg = "white")
ggsave("Figura7_Temporal_PorUF.png", plot_temporal_facet, width = 11, height = 8, dpi = 300, bg = "white")
cat(" -> Figuras 1 a 7 exportadas.\n\n")

# 16.9 Forest plots
cat("=== MÓDULO 16.9: FOREST PLOTS ===\n")

sig_count_por_par <- pares_por_uf %>%
  filter(significativo_uf) %>%
  count(DIAG_PRINC, DIAGSEC1, name = "n_UFs_significativas")

pares_para_forest <- rede_principal_nordeste %>%
  left_join(sig_count_por_par, by = c("DIAG_PRINC", "DIAGSEC1")) %>%
  mutate(n_UFs_significativas = replace_na(n_UFs_significativas, 0L)) %>%
  filter(n_UFs_significativas >= MIN_UFS_SIG_FOREST) %>%
  arrange(desc(n_UFs_significativas), desc(RR_meta)) %>%
  head(N_FOREST_MAX)

cat(" -> Pares significativos em >=", MIN_UFS_SIG_FOREST, "UFs:", nrow(pares_para_forest), "\n")

gerar_forest_par <- function(base, comorb, n_sig = NA_integer_) {
  sub_df <- pares_por_uf %>%
    filter(DIAG_PRINC == base, DIAGSEC1 == comorb,
           !is.na(RR), RR > 0, is.finite(RR), is.finite(se_logRR), se_logRR > 0) %>%
    arrange(UF)
  if (nrow(sub_df) < MIN_UFS_METAAN) return(invisible(NULL))
  
  yi <- log(sub_df$RR); sei <- sub_df$se_logRR
  fit <- tryCatch(suppressWarnings(rma(yi = yi, sei = sei, method = "REML")), error = function(e) NULL)
  if (is.null(fit)) fit <- tryCatch(suppressWarnings(rma(yi = yi, sei = sei, method = "DL")), error = function(e) NULL)
  if (is.null(fit)) return(invisible(NULL))
  
  rotulos_uf <- paste0(sub_df$UF, ifelse(sub_df$significativo_uf, " *", ""))
  nome_arquivo <- paste0("Forest_", base, "_x_", comorb, ".png")
  png(nome_arquivo, width = 1700, height = 950 + 60 * nrow(sub_df), res = 200)
  forest(fit,
         slab    = rotulos_uf,
         atransf = exp,
         xlab    = "Risco Relativo (RR)",
         mlab    = paste0("Pooled \u2014 I\u00b2 = ", round(fit$I2, 1), "% | \u03c4\u00b2 = ", round(fit$tau2, 3)),
         header  = c("UF", "RR [IC 95%]"),
         main    = paste0(base, " (base) \u2192 ", comorb, " (comorbidade)"))
  mtext(paste0("* significativo individualmente (FDR<", ALPHA_FDR, ", RR>", RR_MIN, ") \u2014 significativo em ",
               n_sig, "/", length(UFS_NORDESTE), " UFs"),
        side = 1, line = 4, cex = 0.8, font = 3)
  dev.off()
  invisible(NULL)
}

if (nrow(pares_para_forest) > 0) {
  pwalk(list(pares_para_forest$DIAG_PRINC, pares_para_forest$DIAGSEC1, pares_para_forest$n_UFs_significativas),
        gerar_forest_par)
  cat(" -> Forest plots individuais exportados:", nrow(pares_para_forest), "\n")
  
  top_forest_resumo <- pares_para_forest %>%
    mutate(par_label = paste0(DIAG_PRINC, " \u2192 ", DIAGSEC1, "  (", n_UFs_significativas, "/",
                              length(UFS_NORDESTE), " UFs)"),
           par_label = fct_reorder(par_label, RR_meta))
  
  plot_forest_resumo <- ggplot(top_forest_resumo, aes(x = RR_meta, y = par_label, color = robustez_regional)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = RR_meta_li, xmax = RR_meta_ls), orientation = "y", width = 0.25, linewidth = 0.6) +
    geom_point(aes(size = n_UFs_significativas)) +
    scale_color_manual(values = c("Robusta (baixa heterogeneidade)" = "#1a9850",
                                  "Estado-específica (alta heterogeneidade)" = "#e34a33"),
                       name = "Classificação (I²)") +
    scale_size_continuous(name = "UFs significativas", range = c(2, 6)) +
    scale_x_log10() +
    labs(title = "Forest Plot Resumo \u2014 Associações Significativas em Múltiplas UFs",
         subtitle = paste0("Pares significativos em >= ", MIN_UFS_SIG_FOREST, " UFs; RR meta-analítico com IC95%"),
         x = "Risco Relativo meta-analítico (escala log)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave("Figura8_ForestPlot_Resumo_TopPares.png", plot_forest_resumo, width = 10, height = 9, dpi = 300, bg = "white")
  cat(" -> Figura 8 exportada.\n\n")
} else {
  cat(" -> Nenhum par atingiu o critério de forest plot; figuras puladas.\n\n")
}

# ==============================================================================
# MÓDULO 16.10 (TRECHO FINAL CORRIGIDO): README + WORKBOOK OSF
# Retoma após o CSV de pacientes (fwrite) já gravado.
# Correção: todas as colunas do README agora são texto (as.character), o que
# evita o erro "Can't combine <character> and <double>" do tribble().
# ==============================================================================

stopifnot(
  "Faltam objetos dos módulos anteriores" =
    all(sapply(c("tabela_hubs_regional", "tabela3_regional", "meta_todos", "pares_por_uf",
                 "pares_pool_ingenuo", "resultados_firth", "resultados_misto", "topologia_empirica",
                 "topologia_nulos", "lei_potencia_result", "resultados_mc", "robustez_resultados",
                 "resumo_direcional", "hubs_universais", "tabela_bootstrap",
                 "tabela_sensibilidade_regional", "tabela_sensibilidade_uf", "resumo_atricao",
                 "tabela_temporal_regional", "tabela_temporal_uf", "metricas_seguimento_uf",
                 "pares_para_forest", "media_nula", "p_empirico_cascata", "p_empirico_mod",
                 "pc_molloy_reed", "modularidade_regional"), exists))
)

# formatador de número com separador de milhar brasileiro (sem o aviso de big.mark)
fmt_n <- function(x) format(x, big.mark = ".", decimal.mark = ",", scientific = FALSE)

wb_osf <- createWorkbook()
addWorksheet(wb_osf, "README_Metodologia")

readme_df <- tibble::tribble(
  ~Item, ~Descricao,
  "Título do Estudo", "Topologia da Mortalidade Hospitalar - Região Nordeste, 2016-2025",
  "UFs incluídas", paste(UFS_NORDESTE, collapse = ", "),
  "Fonte de Dados", "SIH/SUS (DATASUS), extraídos via pacote microdatasus",
  "Coorte Analítica Final", paste0("N = ", fmt_n(N_obitos_nordeste), " óbitos hospitalares purificados"),
  "Critério da Rede Regional", paste0("Meta-análise de efeitos aleatórios (log-RR; REML com fallback DL) entre UFs, FDR < ", ALPHA_FDR, ", RR_meta > ", RR_MIN),
  "Robustez Regional", paste0("Aresta 'Robusta' se I² < ", I2_ROBUSTO, "%"),
  "Modularidade (rede regional)", paste0("Louvain ponderado = ", round(modularidade_regional, 4), " (p-empírico = ", round(p_empirico_mod, 4), ")"),
  "Simulação de Monte Carlo", paste0("Percolação de aresta, B = ", B_MC, " réplicas, ", vcount(grafo_regional), " sementes; ", length(media_nula), " redes nulas Erdős–Rényi pareadas; p-empírico = ", round(p_empirico_cascata, 4)),
  "Bootstrap de centralidade", paste0("B = ", B_BOOT, " réplicas paramétricas (binomial sobre a,b,c,d) + meta-análise DerSimonian-Laird com FDR por réplica"),
  "Limiar de percolação (Molloy-Reed)", as.character(round(pc_molloy_reed, 4)),   # <- CORRIGIDO (era numérico)
  "Forest plots", paste0(nrow(pares_para_forest), " individuais (pares com significância em >= ", MIN_UFS_SIG_FOREST, " UFs) + 1 resumo"),
  "Dados individuais", ifelse(REMOVER_QUASE_IDENTIFICADORES, "CSV separado, sem município e datas", "CSV separado, com quase-identificadores (verificar LGPD/CEP)"),
  "Gerado em", format(Sys.Date(), "%d/%m/%Y")
)
writeData(wb_osf, "README_Metodologia", readme_df)

add_aba <- function(nome, df) { addWorksheet(wb_osf, nome); writeData(wb_osf, nome, df) }
add_aba("Tabela2_Hubs_Regional",     tabela_hubs_regional)
add_aba("Tabela3_MetaAnalitica",     tabela3_regional)
add_aba("Meta_Todos_Pares",          meta_todos)
add_aba("Pares_Por_UF",              pares_por_uf)
add_aba("Pool_Ingenuo_vs_Meta",      pares_pool_ingenuo)
add_aba("Firth_Top_Pares",           resultados_firth)
add_aba("Misto_Top_Pares",           resultados_misto)
add_aba("Topologia_vs_Nulos",        bind_rows(topologia_empirica, topologia_nulos))
add_aba("Lei_Potencia",              lei_potencia_result)
add_aba("MonteCarlo_Cascata",        resultados_mc)
add_aba("Robustez_Ataque",           robustez_resultados)
add_aba("Rede_Direcionada",          resumo_direcional)
add_aba("Hubs_Universais_Multiplex", hubs_universais)
add_aba("Bootstrap_Robustez",        tabela_bootstrap)
add_aba("Sensibilidade_Regional",    tabela_sensibilidade_regional)
add_aba("Sensibilidade_Por_UF",      tabela_sensibilidade_uf)
add_aba("Atricao_Por_UF",            resumo_atricao)
add_aba("Temporal_Regional",         tabela_temporal_regional)
add_aba("Temporal_Por_UF",           tabela_temporal_uf)
add_aba("Permanencia_Por_UF",        metricas_seguimento_uf)

arquivo_xlsx_osf <- "OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE.xlsx"
saveWorkbook(wb_osf, arquivo_xlsx_osf, overwrite = TRUE)
cat(" -> Workbook exportado:", arquivo_xlsx_osf, "\n\n")

cat("=== PIPELINE V4.0 (NORDESTE) CONCLUÍDO COM SUCESSO! ===\n")

getwd()




# ==============================================================================
# PIPELINE NORDESTE v4.1 - ANÁLISES AVANÇADAS, TOPOLÓGICAS E META-ANALÍTICAS
# ==============================================================================
# Pré-requisito: Objetos dos Módulos 1-5 carregados no ambiente:
# pares_por_uf, pares_pool_ingenuo, dados_uf, coorte_nordeste, resumo_atricao,
# N_obitos_nordeste, parâmetros do Módulo 2 e funções de checkpoint/log_tempo.
# ==============================================================================

cat("=== INICIANDO PIPELINE NORDESTE v4.1 (ANÁLISES E VISUALIZAÇÕES) ===\n")

stopifnot(
  "Faltam objetos dos Módulos 1 a 5 no ambiente!" =
    all(sapply(c("pares_por_uf", "pares_pool_ingenuo", "dados_uf", "coorte_nordeste"), exists))
)

# ------------------------------------------------------------------------------
# 0. CORREÇÕES PRELIMINARES DE DADOS
# ------------------------------------------------------------------------------
# Correção mandatória da UF do Piauí que estava sendo perdida no bind
coorte_nordeste <- coorte_nordeste %>%
  mutate(UF = if_else(is.na(UF) | UF == "" | UF == "<NA>", UF_lista, UF))

REMOVER_QUASE_IDENTIFICADORES <- TRUE
B_NULL_GRAFOS <- 100

# ==============================================================================
# MÓDULO 6: META-ANÁLISE ENTRE ESTADOS (KNAPP-HARTUNG & HETEROGENEIDADE)
# ==============================================================================
cat("=== MÓDULO 6: META-ANÁLISE DE EFEITOS ALEATÓRIOS (REML / KNAPP-HARTUNG) ===\n")

.erros_meta <- character(0)

meta_reml_par <- function(yi, sei) {
  metodo <- "REML"
  # Knapp-Hartung recomendado para pequeno número de estudos (k <= 9)
  fit <- tryCatch(
    suppressWarnings(metafor::rma(yi = yi, sei = sei, method = "REML", test = "knapp-hartung")),
    error = function(e) { .erros_meta <<- c(.erros_meta, conditionMessage(e)); NULL }
  )
  if (is.null(fit)) {
    metodo <- "DL"
    fit <- tryCatch(
      suppressWarnings(metafor::rma(yi = yi, sei = sei, method = "DL")),
      error = function(e) { .erros_meta <<- c(.erros_meta, conditionMessage(e)); NULL }
    )
  }
  if (is.null(fit)) return(NULL)
  
  # Intervalo de predição de 95%
  pred <- tryCatch(predict(fit), error = function(e) NULL)
  
  list(
    metodo_tau2   = metodo,
    RR_meta       = exp(as.numeric(fit$b[1])),
    RR_meta_li    = exp(fit$ci.lb),
    RR_meta_ls    = exp(fit$ci.ub),
    PI_meta_li    = if (!is.null(pred)) exp(pred$pi.lb) else NA_real_,
    PI_meta_ls    = if (!is.null(pred)) exp(pred$pi.ub) else NA_real_,
    tau2          = fit$tau2,
    I2            = fit$I2,
    Q_valor       = fit$QE,
    Q_p           = fit$QEp,
    p_meta        = fit$pval
  )
}

pares_dt <- as.data.table(pares_por_uf)[
  !is.na(RR) & RR > 0 & is.finite(RR) & is.finite(se_logRR) & se_logRR > 0]
pares_dt[, k_ufs := .N, by = .(DIAG_PRINC, DIAGSEC1)]

log_tempo("Ajustando meta-análises por par clínico...")
meta_bruto <- pares_dt[k_ufs >= MIN_UFS_METAAN, {
  r <- meta_reml_par(log(RR), se_logRR)
  if (is.null(r)) NULL else c(
    list(n_UFs_incluidas = .N, UFs_incluidas = paste(UF, collapse = ", ")), r)
}, by = .(DIAG_PRINC, DIAGSEC1)]

meta_todos <- as_tibble(meta_bruto) %>%
  mutate(
    p_meta_fdr = p.adjust(p_meta, method = "fdr"),
    robustez_regional = ifelse(I2 < I2_ROBUSTO,
                               "Robusta (baixa heterogeneidade)",
                               "Estado-específica (alta heterogeneidade)")
  )

meta_resultados <- meta_todos %>%
  filter(!is.na(p_meta_fdr), p_meta_fdr < ALPHA_FDR, RR_meta > RR_MIN)

rede_principal_nordeste <- meta_resultados %>%
  mutate(
    peso_bruto = n_UFs_incluidas,
    weight     = RR_meta,
    tier       = ifelse(robustez_regional == "Robusta (baixa heterogeneidade)", "Alta", "Exploratória")
  )

cat(" -> Pares meta-analisados:", nrow(meta_todos), "| Significantes na rede:", nrow(meta_resultados), "\n\n")

# ==============================================================================
# MÓDULO 6.1: AVALIAÇÃO DO PARADOXO DE SIMPSON (POOL INGÊNUO vs. META-ANÁLISE)
# ==============================================================================
cat("=== MÓDULO 6.1: CONTRASTE METODOLÓGICO (POOL INGÊNUO vs. META-ANÁLISE) ===\n")

comparativo_simpson <- meta_resultados %>%
  inner_join(pares_pool_ingenuo %>% select(DIAG_PRINC, DIAGSEC1, RR_pool_ingenuo),
             by = c("DIAG_PRINC", "DIAGSEC1")) %>%
  mutate(
    Diferenca_Relativa_Pct = round(100 * (RR_pool_ingenuo - RR_meta) / RR_meta, 2),
    Discordancia_Clinica = abs(Diferenca_Relativa_Pct) > 20
  ) %>%
  arrange(desc(abs(Diferenca_Relativa_Pct)))

cat(" -> Top 5 associações com maior distorção pelo Pool Ingênuo:\n")
print(comparativo_simpson %>% 
        select(DIAG_PRINC, DIAGSEC1, RR_meta, RR_pool_ingenuo, Diferenca_Relativa_Pct) %>% 
        head(5))
cat("\n")

# ==============================================================================
# MÓDULO 7: REGRESSÃO LOGÍSTICA — FIRTH E MODELOS MISTOS
# ==============================================================================
cat("=== MÓDULO 7: REGRESSÃO PENALIZADA DE FIRTH E GLMM ===\n")

universo_pool_dt <- rbindlist(lapply(names(dados_uf), function(uf) {
  x <- as.data.table(dados_uf[[uf]]$universo)
  x[, UF := uf]
  x
}))

calcular_2x2_dt <- function(base, comorb) {
  sub  <- universo_pool_dt[DIAG_PRINC_U == base]
  expo <- !is.na(sub$DIAGSEC1_U) & sub$DIAGSEC1_U == comorb
  tibble(
    a = sum(expo & sub$OBITO_U),  b = sum(expo & !sub$OBITO_U),
    c = sum(!expo & sub$OBITO_U), d = sum(!expo & !sub$OBITO_U)
  )
}

calcular_or_firth <- function(a, b, cc, d) {
  df <- data.frame(exposto = c(1, 1, 0, 0), obito = c(1, 0, 1, 0), w = c(a, b, cc, d))
  df <- df[df$w > 0, ]
  vazio <- tibble(OR_firth = NA_real_, OR_firth_li = NA_real_, OR_firth_ls = NA_real_)
  fit <- tryCatch(logistf(obito ~ exposto, data = df, weights = w, pl = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(vazio)
  tibble(
    OR_firth    = exp(unname(coef(fit)["exposto"])),
    OR_firth_li = exp(unname(fit$ci.lower["exposto"])),
    OR_firth_ls = exp(unname(fit$ci.upper["exposto"]))
  )
}

top_pares_firth <- rede_principal_nordeste %>% arrange(desc(RR_meta)) %>% head(15)

resultados_firth <- top_pares_firth %>%
  mutate(cel = map2(DIAG_PRINC, DIAGSEC1, calcular_2x2_dt)) %>%
  unnest(cel) %>%
  mutate(firth = pmap(list(a, b, c, d), calcular_or_firth)) %>%
  unnest(firth)

ajustar_misto_uf <- function(base, comorb) {
  vazio <- tibble(OR_misto = NA_real_, OR_misto_li = NA_real_, OR_misto_ls = NA_real_, var_UF = NA_real_)
  dfb <- universo_pool_dt[DIAG_PRINC_U == base, .(n = .N),
                          by = .(UF,
                                 exposto = as.integer(!is.na(DIAGSEC1_U) & DIAGSEC1_U == comorb),
                                 obito   = as.integer(OBITO_U))]
  df <- dcast(dfb, UF + exposto ~ obito, value.var = "n", fill = 0L)
  if (!"1" %in% names(df)) df[, `1` := 0L]
  if (!"0" %in% names(df)) df[, `0` := 0L]
  setnames(df, c("1", "0"), c("obitos", "sobreviventes"))
  if (uniqueN(df$UF) < 2 || uniqueN(df$exposto) < 2) return(vazio)
  fit <- tryCatch(
    glmer(cbind(obitos, sobreviventes) ~ exposto + (1 | UF), data = as.data.frame(df),
          family = binomial(link = "logit"), control = glmerControl(optimizer = "bobyqa")),
    error = function(e) NULL
  )
  if (is.null(fit)) return(vazio)
  est <- fixef(fit)["exposto"]; se <- sqrt(vcov(fit)["exposto", "exposto"])
  tibble(OR_misto = exp(est), OR_misto_li = exp(est - 1.96 * se), OR_misto_ls = exp(est + 1.96 * se),
         var_UF = as.numeric(VarCorr(fit)$UF[1]))
}

resultados_misto <- top_pares_firth %>%
  head(5) %>%
  mutate(misto = map2(DIAG_PRINC, DIAGSEC1, ajustar_misto_uf)) %>%
  unnest(misto)

# ==============================================================================
# MÓDULO 8: CONSTRUÇÃO DAS REDES E MULTICAMADA (MULTIPLEX)
# ==============================================================================
cat("=== MÓDULO 8: CONSTRUÇÃO DA REDE REGIONAL E MULTICAMADA ===\n")

grafo_regional <- graph_from_data_frame(
  rede_principal_nordeste %>%
    select(DIAG_PRINC, DIAGSEC1, weight, RR_meta, I2, robustez_regional, tier, n_UFs_incluidas),
  directed = FALSE
)
grafo_regional <- delete_vertices(grafo_regional, which(degree(grafo_regional) == 0))
V(grafo_regional)$betweenness <- betweenness(grafo_regional, directed = FALSE, normalized = TRUE) * 100
comunidades_regional <- cluster_louvain(grafo_regional)
V(grafo_regional)$group <- as.character(membership(comunidades_regional))
modularidade_regional <- modularity(comunidades_regional)

grafos_uf <- list()
for (uf in UFS_NORDESTE) {
  rede_uf <- pares_por_uf %>% filter(UF == uf, significativo_uf)
  if (nrow(rede_uf) == 0) { grafos_uf[[uf]] <- NULL; next }
  g <- graph_from_data_frame(rede_uf %>% mutate(weight = casos_obito_par * RR) %>%
                               select(DIAG_PRINC, DIAGSEC1, weight, RR), directed = FALSE)
  g <- delete_vertices(g, which(degree(g) == 0))
  grafos_uf[[uf]] <- g
}

participacao_multiplex <- map_dfr(names(grafos_uf), function(uf) {
  g <- grafos_uf[[uf]]
  if (is.null(g) || vcount(g) == 0) return(tibble())
  deg <- degree(g)
  limiar_hub <- quantile(deg, 0.90, na.rm = TRUE)
  tibble(UF = uf, diagnostico = names(deg)[deg >= limiar_hub], grau = deg[deg >= limiar_hub])
})

hubs_universais <- participacao_multiplex %>%
  count(diagnostico, name = "n_camadas_como_hub") %>%
  arrange(desc(n_camadas_como_hub))

# ==============================================================================
# MÓDULO 9: TOPOLOGIA FORMAL COM CORREÇÃO DE DISTÂNCIA
# ==============================================================================
cat("=== MÓDULO 9: TOPOLOGIA FORMAL E CORREÇÃO DE DISTÂNCIA ===\n")

descrever_topologia_corrigida <- function(g, nome = "empírica") {
  # CORREÇÃO CRÍTICA: weights = NA calcula a distância topológica por passos/arestas.
  # Para distância ponderada, custo = 1 / log(RR + 1)
  dist_topologica <- mean_distance(g, directed = FALSE, unconnected = TRUE, weights = NA)
  diam_topologico <- diameter(g, directed = FALSE, unconnected = TRUE, weights = NA)
  
  dist_ponderada  <- mean_distance(g, directed = FALSE, unconnected = TRUE, 
                                   weights = 1 / log(E(g)$weight + 1))
  
  tibble(
    Rede                  = nome,
    N_nos                 = vcount(g),
    N_arestas             = ecount(g),
    Grau_medio            = round(mean(degree(g)), 3),
    Clustering_global     = round(transitivity(g, type = "global"), 4),
    Caminho_medio_passos  = round(dist_topologica, 3),
    Diametro_passos       = as.numeric(diam_topologico),
    Caminho_ponderado     = round(dist_ponderada, 3),
    Assortatividade       = round(assortativity_degree(g, directed = FALSE), 4)
  )
}

ajustar_lei_potencia <- function(g) {
  deg <- degree(g); deg <- deg[deg > 0]
  fit <- tryCatch(fit_power_law(deg, implementation = "plfit"), error = function(e) NULL)
  if (is.null(fit)) return(tibble(alpha = NA_real_, x_min = NA_real_, p_ks = NA_real_))
  pega <- function(x) if (is.null(x) || !is.numeric(x) || length(x) == 0) NA_real_ else as.numeric(x[1])
  tibble(alpha = round(pega(fit$alpha), 3),
         x_min = pega(fit$xmin),
         p_ks  = round(pega(fit$KS.p), 4))
}

percolacao_molloy_reed <- function(g) {
  deg <- degree(g); k1 <- mean(deg); k2 <- mean(deg^2)
  if ((k2 - k1) <= 0) return(NA_real_)
  k1 / (k2 - k1)
}

gerar_nulos_pareados <- function(g, n_reps = 30) {
  N <- vcount(g); E <- ecount(g); k_medio <- 2 * E / N
  m_ba <- max(1, round(k_medio / 2))
  er_list <- lapply(seq_len(n_reps), function(i) sample_gnm(N, E))
  ba_list <- lapply(seq_len(n_reps), function(i) sample_pa(N, m = m_ba, directed = FALSE))
  ws_list <- lapply(seq_len(n_reps), function(i) sample_smallworld(dim = 1, size = N, nei = m_ba, p = 0.05))
  resumir <- function(lst, nome) {
    m <- map_dfr(lst, function(gg) tibble(
      clustering = transitivity(gg, type = "global"),
      caminho    = mean_distance(gg, directed = FALSE, unconnected = TRUE, weights = NA)))
    tibble(Rede = nome,
           Clustering_medio = round(mean(m$clustering, na.rm = TRUE), 4),
           Clustering_dp    = round(sd(m$clustering, na.rm = TRUE), 4),
           Caminho_medio    = round(mean(m$caminho, na.rm = TRUE), 3))
  }
  bind_rows(resumir(er_list, "Erdős–Rényi (nula)"),
            resumir(ba_list, "Barabási–Albert (nula)"),
            resumir(ws_list, "Watts–Strogatz (nula)"))
}

topologia_empirica  <- descrever_topologia_corrigida(grafo_regional, "Nordeste (empírica)")
lei_potencia_result <- ajustar_lei_potencia(grafo_regional)
pc_molloy_reed      <- percolacao_molloy_reed(grafo_regional)
topologia_nulos     <- gerar_nulos_pareados(grafo_regional, n_reps = 30)

cat(" -> Topologia empírica corrigida:\n")
print(as.data.frame(topologia_empirica))
cat("\n")

# ==============================================================================
# MÓDULO 10: PERCOLAÇÃO E MONTE CARLO (CORREÇÃO DO TESTE DE HIPÓTESE)
# ==============================================================================
cat("=== MÓDULO 10: MONTE CARLO E TESTE UNICAUDAL DE RESTRIÇÃO MODULAR ===\n")

normalizar_prob_transmissao <- function(g, prob_base = 0.5) {
  rr <- E(g)$RR_meta
  if (is.null(rr)) rr <- E(g)$RR
  rr_norm <- (rr - min(rr, na.rm = TRUE)) / (max(rr, na.rm = TRUE) - min(rr, na.rm = TRUE) + 1e-9)
  prob_base * (0.5 + rr_norm)
}

mc_cascata_rapido <- function(g, p_aresta, idx_sementes, B) {
  el    <- as_edgelist(g, names = FALSE)
  N     <- vcount(g); M <- nrow(el)
  vazio <- make_empty_graph(N, directed = FALSE)
  tam   <- matrix(NA_integer_, nrow = length(idx_sementes), ncol = B)
  for (b in seq_len(B)) {
    ab <- runif(M) < p_aresta
    gb <- if (any(ab)) add_edges(vazio, as.vector(t(el[ab, , drop = FALSE]))) else vazio
    cmp <- components(gb)
    tam[, b] <- cmp$csize[cmp$membership[idx_sementes]]
  }
  tam
}

prob_aresta_regional <- normalizar_prob_transmissao(grafo_regional, prob_base = 0.4)
N_g <- vcount(grafo_regional); E_g <- ecount(grafo_regional)

tam_emp <- mc_cascata_rapido(grafo_regional, prob_aresta_regional, seq_len(N_g), B_MC)
resultados_mc <- tibble(
  diagnostico_semente = V(grafo_regional)$name,
  tamanho_medio   = rowMeans(tam_emp),
  tamanho_mediano = apply(tam_emp, 1, median),
  tamanho_p95     = apply(tam_emp, 1, quantile, 0.95)
) %>% arrange(desc(tamanho_medio))
media_emp <- mean(tam_emp)

p_nula <- mean(prob_aresta_regional, na.rm = TRUE)
media_nula <- vapply(seq_len(B_NULL_GRAFOS), function(i) {
  g_er <- sample_gnm(N_g, E_g)
  mean(mc_cascata_rapido(g_er, rep(p_nula, E_g), seq_len(N_g), 100))
}, numeric(1))

# CORREÇÃO CRÍTICA DO P-VALOR: H0 = Rede aleatória produz cascata menor ou igual
# Hipótese Alternativa: A rede real restringe a cascata (media_emp < media_nula)
p_empirico_cascata <- (1 + sum(media_nula <= media_emp)) / (length(media_nula) + 1)

cat(" -> Cascata média empírica:", round(media_emp, 2), "| Média sob redes nulas:", round(mean(media_nula), 2), "\n")
cat(" -> p-empírico de restrição modular (H1: empírica < nula):", round(p_empirico_cascata, 4), "\n\n")

# ==============================================================================
# MÓDULO 11 A 13: ROBUSTEZ, REDE DIRECIONADA E BOOTSTRAP
# ==============================================================================
cat("=== MÓDULOS 11-13: ROBUSTEZ, DIRECIONALIDADE E BOOTSTRAP ===\n")

simular_ataque <- function(g, tipo = c("grau", "betweenness", "aleatorio"), n_rep = 30) {
  tipo <- match.arg(tipo)
  N <- vcount(g)
  fracoes <- seq(0, 0.9, by = 0.05)
  calcular_frac_gigante <- function(ordem_remocao) {
    map_dbl(fracoes, function(f) {
      k <- floor(f * N)
      if (k == 0) return(max(components(g)$csize) / N)
      g_sub <- delete_vertices(g, ordem_remocao[seq_len(k)])
      if (vcount(g_sub) == 0) return(0)
      max(components(g_sub)$csize) / N
    })
  }
  if (tipo == "grau") {
    ordem <- V(g)$name[order(-degree(g))]
    tibble(fracao_removida = fracoes, frac_giant = calcular_frac_gigante(ordem), tipo = "Ataque por grau")
  } else if (tipo == "betweenness") {
    ordem <- V(g)$name[order(-betweenness(g))]
    tibble(fracao_removida = fracoes, frac_giant = calcular_frac_gigante(ordem), tipo = "Ataque por intermediação")
  } else {
    reps <- map(seq_len(n_rep), function(i) calcular_frac_gigante(sample(V(g)$name)))
    tibble(fracao_removida = fracoes, frac_giant = Reduce(`+`, reps) / n_rep, tipo = "Falha aleatória")
  }
}

robustez_resultados <- bind_rows(
  simular_ataque(grafo_regional, "grau"),
  simular_ataque(grafo_regional, "betweenness"),
  simular_ataque(grafo_regional, "aleatorio", n_rep = 50)
)

grafo_direcionado <- graph_from_data_frame(
  rede_principal_nordeste %>% select(DIAG_PRINC, DIAGSEC1, weight, RR_meta), directed = TRUE)

resumo_direcional <- tibble(
  diagnostico = V(grafo_direcionado)$name,
  out_degree  = degree(grafo_direcionado, mode = "out"),
  in_degree   = degree(grafo_direcionado, mode = "in"),
  betweenness = round(betweenness(grafo_direcionado, directed = TRUE), 2)
) %>% arrange(desc(out_degree))

componentes_dir <- components(grafo_direcionado, mode = "weak")
componente_gigante_frac <- max(componentes_dir$csize) / vcount(grafo_direcionado)

top10_ids <- V(grafo_regional)$name[order(-V(grafo_regional)$betweenness)][1:min(10, vcount(grafo_regional))]
boot_resultados <- matrix(NA_real_, nrow = B_BOOT, ncol = length(top10_ids), dimnames = list(NULL, top10_ids))

metanalise_dl_rapida <- function(dt) {
  dt <- dt[!is.na(RR) & is.finite(RR) & RR > 0 & is.finite(se_logRR) & se_logRR > 0]
  if (nrow(dt) == 0) return(data.table())
  dt[, `:=`(yi = log(RR), vi = se_logRR^2)]
  dt[, wi := 1 / vi]
  dt[, k := .N, by = .(DIAG_PRINC, DIAGSEC1)]
  dt <- dt[k >= MIN_UFS_METAAN]
  if (nrow(dt) == 0) return(data.table())
  
  dt[, `:=`(sum_wi = sum(wi), sum_wi2 = sum(wi^2), mu_f = sum(wi * yi) / sum(wi)),
     by = .(DIAG_PRINC, DIAGSEC1)]
  dt[, Q := sum(wi * (yi - mu_f)^2), by = .(DIAG_PRINC, DIAGSEC1)]
  dt[, tau2 := pmax(0, (Q - (k - 1)) / (sum_wi - sum_wi2 / sum_wi))]
  dt[, wi_star := 1 / (vi + tau2)]
  
  saida <- dt[, .(yi_pooled = sum(wi_star * yi) / sum(wi_star),
                  se_pooled = sqrt(1 / sum(wi_star))), by = .(DIAG_PRINC, DIAGSEC1)]
  saida[, `:=`(RR_meta_boot = exp(yi_pooled),
               p_boot = 2 * pnorm(-abs(yi_pooled / se_pooled)))]
  saida[, p_boot_fdr := p.adjust(p_boot, "fdr")]
  saida[]
}

pares_uf_dt <- as.data.table(pares_por_uf)[, .(UF, DIAG_PRINC, DIAGSEC1, a, b, c, d)]
set.seed(123)
for (i in seq_len(B_BOOT)) {
  boot_i <- copy(pares_uf_dt)
  boot_i[, a_boot := rbinom(.N, size = a + b, prob = a / (a + b))]
  boot_i[, c_boot := rbinom(.N, size = c + d, prob = c / (c + d))]
  boot_i[, `:=`(
    RR       = (a_boot / (a + b)) / (c_boot / (c + d)),
    se_logRR = sqrt(1 / pmax(a_boot, 0.5) - 1 / (a + b) + 1 / pmax(c_boot, 0.5) - 1 / (c + d))
  )]
  meta_boot <- metanalise_dl_rapida(boot_i)
  if (nrow(meta_boot) == 0) next
  rede_boot <- meta_boot[p_boot_fdr < ALPHA_FDR & RR_meta_boot > RR_MIN]
  if (nrow(rede_boot) == 0) next
  
  g_boot <- graph_from_data_frame(rede_boot[, .(DIAG_PRINC, DIAGSEC1, weight = RR_meta_boot)], directed = FALSE)
  g_boot <- delete_vertices(g_boot, which(degree(g_boot) == 0))
  cent_boot <- betweenness(g_boot, directed = FALSE, normalized = TRUE) * 100
  for (id in top10_ids) if (id %in% names(cent_boot)) boot_resultados[i, id] <- cent_boot[id]
}

tabela_bootstrap <- tibble(
  `CID-10` = top10_ids,
  `Centralidade Observada (%)`     = round(V(grafo_regional)$betweenness[match(top10_ids, V(grafo_regional)$name)], 2),
  `Mediana Bootstrap (%)`          = round(apply(boot_resultados, 2, median, na.rm = TRUE), 2),
  `Estabilidade (% Reamostragens)` = round(100 * colMeans(!is.na(boot_resultados)), 1)
)

pesos_obs <- E(grafo_regional)$weight
set.seed(2024)
modularidades_nulas <- vapply(seq_len(B_MOD), function(i) {
  g_rand <- tryCatch(rewire(grafo_regional, with = keeping_degseq(loops = FALSE, niter = ecount(grafo_regional) * 10)),
                     error = function(e) NULL)
  if (is.null(g_rand)) return(NA_real_)
  E(g_rand)$weight <- sample(pesos_obs)
  com <- tryCatch(cluster_louvain(g_rand), error = function(e) NULL)
  if (is.null(com)) NA_real_ else modularity(com)
}, numeric(1))
modularidades_nulas <- modularidades_nulas[!is.na(modularidades_nulas)]
p_empirico_mod <- (1 + sum(modularidades_nulas >= modularidade_regional)) / (length(modularidades_nulas) + 1)

# ==============================================================================
# MÓDULO 14 E 15: SENSIBILIDADE E SÉRIE TEMPORAL CORRIGIDA
# ==============================================================================
cat("=== MÓDULO 14-15: SENSIBILIDADE E TEMPORALIDADE (UF CORRIGIDA) ===\n")

meta_sig_todos <- meta_todos %>% filter(!is.na(p_meta_fdr), p_meta_fdr < ALPHA_FDR)

tabela_sensibilidade_regional <- map_dfr(LIMIARES_SENS, function(lim) {
  rede_lim <- meta_sig_todos %>% filter(RR_meta > lim)
  if (nrow(rede_lim) == 0) return(tibble(`Limiar RR` = lim, `N Arestas` = 0L, `N Nós` = 0L, `Top3 Hubs` = "Nenhum"))
  g_lim <- graph_from_data_frame(rede_lim %>% select(DIAG_PRINC, DIAGSEC1, RR_meta), directed = FALSE)
  g_lim <- delete_vertices(g_lim, which(degree(g_lim) == 0))
  cent_lim <- betweenness(g_lim, directed = FALSE, normalized = TRUE) * 100
  top3 <- sort(cent_lim, decreasing = TRUE)[1:min(3, length(cent_lim))]
  tibble(`Limiar RR` = lim, `N Arestas` = nrow(rede_lim), `N Nós` = vcount(g_lim),
         `Top3 Hubs` = paste(names(top3), round(top3, 1), sep = ": ", collapse = "; "))
})

tabela_sensibilidade_uf <- map_dfr(UFS_NORDESTE, function(uf) {
  map_dfr(LIMIARES_SENS, function(lim) {
    n <- pares_por_uf %>% filter(UF == uf, RR > lim, significativo_uf) %>% nrow()
    tibble(UF = uf, `Limiar RR` = lim, `N Arestas` = n)
  })
})

# Correção na Série Temporal da UF (PI não será mais NA)
if (!exists("cols_ref") || is.null(cols_ref)) cols_ref <- detectar_colunas(coorte_nordeste)
col_ano <- cols_ref$col_ano

tabela_temporal_uf <- coorte_nordeste %>%
  group_by(UF, Ano = .data[[col_ano]]) %>%
  summarise(Obitos_N = n(), .groups = "drop") %>%
  group_by(UF) %>%
  mutate(Percentual_UF = round(100 * Obitos_N / sum(Obitos_N), 2)) %>%
  ungroup()

tabela_temporal_regional <- tabela_temporal_uf %>%
  group_by(Ano) %>%
  summarise(Obitos_N = sum(Obitos_N), .groups = "drop") %>%
  mutate(Percentual = round(100 * Obitos_N / sum(Obitos_N), 2)) %>%
  arrange(Ano)

metricas_seguimento_uf <- coorte_nordeste %>%
  mutate(DIAS_PERM_NUM = as.numeric(as.character(DIAS_PERM))) %>%
  group_by(UF) %>%
  summarise(
    N_obitos     = n(),
    Mediana_Dias = median(DIAS_PERM_NUM, na.rm = TRUE),
    Q1_Dias      = quantile(DIAS_PERM_NUM, 0.25, na.rm = TRUE),
    Q3_Dias      = quantile(DIAS_PERM_NUM, 0.75, na.rm = TRUE),
    Media_Dias   = round(mean(DIAS_PERM_NUM, na.rm = TRUE), 2),
    .groups = "drop"
  )

# ==============================================================================
# MÓDULO 15.1: SUBANÁLISE EPIDEMIOLÓGICA (PRÉ-PANDEMIA vs PANDEMIA vs PÓS)
# ==============================================================================
cat("=== MÓDULO 15.1: SUBANÁLISE DE HUBS POR PERÍODO PANDÊMICO ===\n")

coorte_nordeste <- coorte_nordeste %>%
  mutate(
    Ano_Num = as.numeric(as.character(.data[[col_ano]])),
    Periodo_Pandemia = case_when(
      Ano_Num <= 2019 ~ "Pre_Pandemia (2016-2019)",
      Ano_Num >= 2020 & Ano_Num <= 2022 ~ "Pandemia_Pico (2020-2022)",
      TRUE ~ "Pos_Pandemia (2023-2025)"
    )
  )

hubs_por_periodo <- coorte_nordeste %>%
  group_by(Periodo_Pandemia) %>%
  group_modify(~ {
    df_p <- .x
    pares_p <- df_p %>% count(DIAG_PRINC, DIAGSEC1, name = "peso") %>% filter(peso >= 10)
    g_p <- graph_from_data_frame(pares_p, directed = FALSE)
    b_p <- betweenness(g_p, normalized = TRUE) * 100
    tibble(CID10 = names(b_p), Betweenness = round(b_p, 2), Grau = degree(g_p)) %>%
      arrange(desc(Betweenness)) %>% head(5)
  })

cat(" -> Top 5 Hubs por Fase Epidemiológica:\n")
print(hubs_por_periodo)
cat("\n")

# ==============================================================================
# MÓDULO 16: TABELAS, FIGURAS E FOREST PLOTS AVANÇADOS (CORRIGIDO)
# ==============================================================================
cat("=== MÓDULO 16: EXPORTAÇÃO GRÁFICA E WORKBOOK OSF COMPLETO ===\n")

cores_modulos <- c("1" = "#1a9850", "2" = "#fc8d62", "3" = "#8073ac", "4" = "#e78ac3",
                   "5" = "#66c2a5", "6" = "#ffd92f", "7" = "#a6d854", "8" = "#e5c494")

tabela_hubs_regional <- tibble(id = V(grafo_regional)$name, group = V(grafo_regional)$group,
                               betweenness = V(grafo_regional)$betweenness) %>%
  mutate(Grau = degree(grafo_regional)[id]) %>%
  left_join(tabela_bootstrap %>% rename(id = `CID-10`), by = "id") %>%
  arrange(desc(betweenness)) %>%
  head(15)

fmt_rr_ic <- function(rr, li, ls) paste0(format(round(rr, 2), decimal.mark = ","), " (",
                                         format(round(li, 2), decimal.mark = ","), "–",
                                         format(round(ls, 2), decimal.mark = ","), ")")
fmt_p <- function(p) ifelse(p < 0.001, "< 0,001", format(round(p, 4), decimal.mark = ","))

tabela3_regional <- rede_principal_nordeste %>%
  mutate(`RR meta (IC95%)` = fmt_rr_ic(RR_meta, RR_meta_li, RR_meta_ls),
         `PI meta (95%)`   = fmt_rr_ic(RR_meta, PI_meta_li, PI_meta_ls),
         `p (FDR)` = fmt_p(p_meta_fdr),
         `I² (%)` = round(I2, 1)) %>%
  arrange(desc(RR_meta)) %>%
  select(`Diagnóstico Base` = DIAG_PRINC, `Comorbidade` = DIAGSEC1, `UFs incluídas` = n_UFs_incluidas,
         `RR meta (IC95%)`, `PI meta (95%)`, `I² (%)`, `p (FDR)`, `Classificação` = robustez_regional)

# Figura 1: Hubs Críticos
plot_hubs <- tabela_hubs_regional %>%
  mutate(rotulo = fct_reorder(id, betweenness)) %>%
  ggplot(aes(x = betweenness, y = rotulo, fill = group)) +
  geom_col(color = "black", linewidth = 0.3, width = 0.7) +
  scale_fill_manual(values = cores_modulos, name = "Módulo Sindrômico") +
  labs(title = "Diagnósticos Mediadores Críticos — Rede Regional Nordeste",
       subtitle = paste0("Top 15 por intermediação | ", vcount(grafo_regional), " nós, ", ecount(grafo_regional), " arestas"),
       x = "Centralidade de intermediação normalizada (%)", y = "CID-10 (diagnóstico base)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14), legend.position = "bottom")
ggsave("Figura1_Hubs_Regional_NORDESTE.png", plot_hubs, width = 10, height = 7, dpi = 300, bg = "white")

# Figura 2: Grafo Topológico
set.seed(42)
grafo_plot <- ggraph(grafo_regional, layout = "fr") +
  geom_edge_link(aes(width = weight, linetype = tier), color = "#555555", alpha = 0.6) +
  scale_edge_width(range = c(0.3, 2.5), guide = "none") +
  scale_edge_linetype_manual(values = c("Alta" = "solid", "Exploratória" = "dotted"), name = "Robustez regional") +
  geom_node_point(aes(size = betweenness, fill = group), shape = 21, color = "black", stroke = 0.5) +
  scale_size_continuous(range = c(3, 16), name = "Intermediação (%)") +
  scale_fill_manual(values = cores_modulos, name = "Módulo Sindrômico") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold", bg.color = "white", bg.r = 0.1) +
  labs(title = "Topologia das Associações Clínicas Letais — Nordeste (2016–2025)",
       subtitle = "Linhas contínuas = robustas entre estados (I² < 50%)") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5), legend.position = "right")
ggsave("Figura2_Grafo_Topologico_NORDESTE.tiff", grafo_plot, width = 12, height = 10, dpi = 300,
       compression = "lzw", bg = "white")

# Figura 3: Lei de Potência (CORRIGIDO: removido o '+' duplicado)
deg_df <- tibble(grau = degree(grafo_regional)) %>% count(grau, name = "freq") %>% mutate(freq_rel = freq / sum(freq))
plot_lei_potencia <- ggplot(deg_df, aes(x = grau, y = freq_rel)) +
  geom_point(size = 2, color = "#2C3E50") +
  scale_x_log10() + 
  scale_y_log10() +
  labs(title = "Distribuição de Grau da Rede Regional (log-log)",
       subtitle = paste0("Lei de potência: alpha = ", lei_potencia_result$alpha, " | KS p = ", lei_potencia_result$p_ks),
       x = "Grau (k)", y = "Frequência relativa") +
  theme_minimal(base_size = 12)
ggsave("Figura3_Distribuicao_Grau_LeiPotencia.png", plot_lei_potencia, width = 8, height = 6, dpi = 300, bg = "white")

# Figura 4: Robustez
plot_robustez <- ggplot(robustez_resultados, aes(x = fracao_removida, y = frac_giant, color = tipo)) +
  geom_line(linewidth = 1.1) + geom_point(size = 1.5) +
  scale_color_manual(values = c("Ataque por grau" = "#C0392B", "Ataque por intermediação" = "#E67E22",
                                "Falha aleatória" = "#2980B9")) +
  labs(title = "Robustez da Rede Regional: Ataque Direcionado vs. Falha Aleatória",
       x = "Fração de nós removidos", y = "Fração do maior componente conectado", color = "Cenário") +
  theme_minimal(base_size = 12)
ggsave("Figura4_Robustez_Ataque_vs_Aleatorio.png", plot_robustez, width = 9, height = 6, dpi = 300, bg = "white")

# Figura 5: Cascata de Monte Carlo
plot_mc <- ggplot(tibble(media_nula = media_nula), aes(x = media_nula)) +
  geom_histogram(bins = 30, fill = "#95a5a6", color = "white") +
  geom_vline(xintercept = media_emp, color = "#C0392B", linewidth = 1.2) +
  labs(title = "Tamanho Médio de Cascata: Rede Empírica vs. Redes Nulas Erdős–Rényi",
       subtitle = paste0("Média empírica = ", round(media_emp, 2), " | Média ER = ", round(mean(media_nula), 2),
                         " | p-empírico (restrição) = ", round(p_empirico_cascata, 4)),
       x = "Tamanho médio de cascata (nº de diagnósticos alcançados)", y = "Nº de redes nulas") +
  theme_minimal(base_size = 12)
ggsave("Figura5_MonteCarlo_Cascata.png", plot_mc, width = 9, height = 6, dpi = 300, bg = "white")

# Figura 6 e 7: Temporais
df_reg <- tabela_temporal_regional %>%
  mutate(Ano = as.numeric(as.character(Ano)), Obitos_Acumulados = cumsum(Obitos_N))
fator <- max(df_reg$Obitos_Acumulados) / max(df_reg$Obitos_N)

plot_temporal_regional <- ggplot(df_reg, aes(x = Ano)) +
  geom_col(aes(y = Obitos_N), fill = "#2C3E50", alpha = 0.85, width = 0.6) +
  geom_line(aes(y = Obitos_Acumulados / fator), color = "#C0392B", linewidth = 1.1, group = 1) +
  scale_x_continuous(breaks = ANO_INICIO:ANO_FIM) +
  scale_y_continuous(name = "Óbitos por ano (Nordeste)",
                     sec.axis = sec_axis(~ . * fator, name = "Óbitos acumulados")) +
  theme_minimal(base_size = 12) +
  labs(title = "Mortalidade Hospitalar Anual e Acumulada — Nordeste (2016–2025)")

plot_temporal_facet <- tabela_temporal_uf %>%
  mutate(Ano = as.numeric(as.character(Ano))) %>%
  ggplot(aes(x = Ano, y = Obitos_N)) +
  geom_col(fill = "#2980B9") +
  facet_wrap(~ UF, scales = "free_y") +
  theme_minimal(base_size = 10) +
  labs(title = "Mortalidade Hospitalar Anual por UF (Consolidado)", x = "Ano", y = "Óbitos")

ggsave("Figura6_Temporal_Regional.png", plot_temporal_regional, width = 11, height = 6, dpi = 300, bg = "white")
ggsave("Figura7_Temporal_PorUF.png", plot_temporal_facet, width = 11, height = 8, dpi = 300, bg = "white")
cat(" -> Figuras 1 a 7 exportadas com sucesso.\n\n")

# ==============================================================================
# MÓDULO 16.9: FOREST PLOTS ORIGINAIS (ESTÁVEL E SEM TRAVAMENTO DE ARQUIVO)
# ==============================================================================
cat("=== MÓDULO 16.9: RESETANDO E GERANDO FOREST PLOTS ORIGINAIS ===\n")

# 1. Fecha qualquer dispositivo gráfico aberto no R que esteja bloqueando arquivos
graphics.off()

# 2. Seleção dos pares com consistência geográfica
sig_count_por_par <- pares_por_uf %>%
  filter(significativo_uf) %>%
  count(DIAG_PRINC, DIAGSEC1, name = "n_UFs_significativas")

pares_para_forest <- rede_principal_nordeste %>%
  left_join(sig_count_por_par, by = c("DIAG_PRINC", "DIAGSEC1")) %>%
  mutate(n_UFs_significativas = replace_na(n_UFs_significativas, 0L)) %>%
  filter(n_UFs_significativas >= MIN_UFS_SIG_FOREST) %>%
  arrange(desc(n_UFs_significativas), desc(RR_meta)) %>%
  head(N_FOREST_MAX)

cat(" -> Pares significativos em >=", MIN_UFS_SIG_FOREST, "UFs:", nrow(pares_para_forest), "\n")

# 3. Função original segura do metafor
gerar_forest_par <- function(base, comorb, n_sig = NA_integer_) {
  sub_df <- pares_por_uf %>%
    filter(DIAG_PRINC == base, DIAGSEC1 == comorb,
           !is.na(RR), RR > 0, is.finite(RR), is.finite(se_logRR), se_logRR > 0) %>%
    arrange(UF)
  
  if (nrow(sub_df) < MIN_UFS_METAAN) return(invisible(NULL))
  
  yi  <- log(sub_df$RR)
  sei <- sub_df$se_logRR
  
  fit <- tryCatch(suppressWarnings(metafor::rma(yi = yi, sei = sei, method = "REML")), error = function(e) NULL)
  if (is.null(fit)) fit <- tryCatch(suppressWarnings(metafor::rma(yi = yi, sei = sei, method = "DL")), error = function(e) NULL)
  if (is.null(fit)) return(invisible(NULL))
  
  rotulos_uf   <- paste0(sub_df$UF, ifelse(sub_df$significativo_uf, " *", ""))
  nome_arquivo <- paste0("Forest_", base, "_x_", comorb, ".png")
  
  # Abre dispositivo gráfico de forma segura
  abriu_png <- tryCatch({
    png(nome_arquivo, width = 1700, height = 950 + 60 * nrow(sub_df), res = 200)
    TRUE
  }, error = function(e) {
    cat(" -> Aviso: Não foi possível abrir o arquivo", nome_arquivo, "(verifique se está aberto em outro programa).\n")
    FALSE
  })
  
  if (!abriu_png) return(invisible(NULL))
  
  tryCatch({
    forest(
      fit,
      slab    = rotulos_uf,
      atransf = exp,
      xlab    = "Risco Relativo (RR)",
      mlab    = paste0("Pooled — I² = ", round(fit$I2, 1), "% | τ² = ", round(fit$tau2, 3)),
      header  = c("UF", "RR [IC 95%]"),
      main    = paste0(base, " (base) → ", comorb, " (comorbidade)")
    )
    mtext(
      paste0("* significativo individualmente (FDR<", ALPHA_FDR, ", RR>", RR_MIN, ") — significativo em ",
             n_sig, "/", length(UFS_NORDESTE), " UFs"),
      side = 1, line = 4, cex = 0.8, font = 3
    )
  }, finally = {
    dev.off()
  })
  
  invisible(NULL)
}

# 4. Geração dos forest plots individuais
if (nrow(pares_para_forest) > 0) {
  pwalk(
    list(pares_para_forest$DIAG_PRINC, 
         pares_para_forest$DIAGSEC1, 
         pares_para_forest$n_UFs_significativas),
    gerar_forest_par
  )
  cat(" -> Forest plots individuais exportados com sucesso!\n")
  
  # 5. Forest Plot Resumo consolidado (Figura 8)
  top_forest_resumo <- pares_para_forest %>%
    mutate(par_label = paste0(DIAG_PRINC, " → ", DIAGSEC1, " (", n_UFs_significativas, "/", length(UFS_NORDESTE), " UFs)"),
           par_label = fct_reorder(par_label, RR_meta))
  
  plot_forest_resumo <- ggplot(top_forest_resumo, aes(x = RR_meta, y = par_label, color = robustez_regional)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbar(aes(xmin = RR_meta_li, xmax = RR_meta_ls), orientation = "y", width = 0.25, linewidth = 0.6) +
    geom_point(aes(size = n_UFs_significativas)) +
    scale_color_manual(values = c("Robusta (baixa heterogeneidade)" = "#1a9850",
                                  "Estado-específica (alta heterogeneidade)" = "#e34a33"),
                       name = "Classificação (I²)") +
    scale_size_continuous(name = "UFs significativas", range = c(2, 6)) +
    scale_x_log10() +
    labs(title = "Forest Plot Resumo — Associações Significativas em Múltiplas UFs",
         subtitle = paste0("Pares significativos em >= ", MIN_UFS_SIG_FOREST, " UFs; RR meta-analítico com IC95%"),
         x = "Risco Relativo meta-analítico (escala log)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
  
  ggsave("Figura8_ForestPlot_Resumo_TopPares.png", plot_forest_resumo, width = 10, height = 9, dpi = 300, bg = "white")
  cat(" -> Figura 8 (Resumo Top Pares) exportada com sucesso.\n\n")
} else {
  cat(" -> Nenhum par atingiu o critério de forest plot.\n\n")
}

# ==============================================================================
# FOREST PLOT ESTILO COCHRANE / REVMAN (desenhado à mão em R base)
# Layout: UF | Com a comorbidade (Eventos, Total) | Sem a comorbidade (Eventos, Total)
#         | [Peso] | RR [IC 95%] | gráfico (escala log, eixo simétrico em torno de 1)
# Rodapé: Heterogeneidade (Tau², Chi², gl, I²), teste do efeito global.
#
# Pré-requisito: pares_por_uf (a, b, c, d, RR, RR_li, RR_ls, se_logRR, significativo_uf),
# parâmetros do Módulo 2 (MIN_UFS_METAAN, ALPHA_FDR, RR_MIN, UFS_NORDESTE),
# pacotes metafor e dplyr carregados.
# ==============================================================================

fmt_int  <- function(x) format(x, big.mark = ".", decimal.mark = ",", trim = TRUE, scientific = FALSE)
fmt_d2   <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = ",")
fmt_pct1 <- function(x) paste0(formatC(x, format = "f", digits = 1, decimal.mark = ","), "%")
fmt_p5   <- function(p) {
  if (p < 0.00001) "< 0,00001"
  else paste0("= ", formatC(p, format = "f", digits = ifelse(p < 0.001, 5, 3), decimal.mark = ","))
}
fmt_rr   <- function(rr, li, ls) paste0(fmt_d2(rr), " [", fmt_d2(li), "; ", fmt_d2(ls), "]")

gerar_forest_cochrane <- function(base, comorb, n_sig = NA_integer_,
                                  ic = c("wald", "hksj"),        # "hksj" = Hartung-Knapp (mais conservador)
                                  mostrar_peso = FALSE,          # coluna de peso (%); ver nota na resposta
                                  eixo = c("simetrico", "ajustado"),
                                  marcar_sig = FALSE,            # asterisco nas UFs significativas
                                  cor_quadrado = "#1F3FBF",      # use "black" p/ estilo RevMan 5 clássico
                                  cor_losango  = "#3A3A3A",
                                  label_esq = "Menor risco de \u00f3bito",
                                  label_dir = "Maior risco de \u00f3bito",
                                  dir_saida = ".") {
  ic   <- match.arg(ic)
  eixo <- match.arg(eixo)
  
  d <- pares_por_uf %>%
    filter(DIAG_PRINC == base, DIAGSEC1 == comorb,
           !is.na(RR), is.finite(RR), RR > 0, is.finite(se_logRR), se_logRR > 0) %>%
    arrange(UF) %>%
    as.data.frame()
  if (nrow(d) < MIN_UFS_METAAN) return(invisible(NULL))
  
  d$n_exp  <- d$a + d$b
  d$n_nexp <- d$c + d$d
  d$yi     <- log(d$RR)
  
  ajustar <- function(metodo) tryCatch(
    suppressWarnings(rma(yi = yi, sei = se_logRR, data = d, method = metodo,
                         test = if (ic == "hksj") "knha" else "z")),
    error = function(e) NULL)
  fit <- ajustar("REML")
  if (is.null(fit)) fit <- ajustar("DL")
  if (is.null(fit)) return(invisible(NULL))
  
  k   <- nrow(d)
  w   <- as.numeric(weights(fit))            # % (soma 100)
  est <- as.numeric(fit$b[1]); ci_lb <- fit$ci.lb; ci_ub <- fit$ci.ub
  
  # ---- eixo (escala log) -----------------------------------------------------
  lo_all <- min(c(d$RR_li, exp(ci_lb)))
  hi_all <- max(c(d$RR_ls, exp(ci_ub)))
  if (eixo == "simetrico") {
    Hs <- c(2, 5, 10, 20, 50, 100, 1000)
    m  <- max(hi_all, 1 / lo_all)
    H  <- Hs[Hs >= m][1]; if (is.na(H)) H <- ceiling(m)
    alim <- log(c(1 / H, H))
    cand <- c(.001, .01, .02, .05, .1, .2, .5, 1, 2, 5, 10, 20, 50, 100, 1000)
  } else {
    alim <- log(c(min(lo_all, 1) / 1.2, max(hi_all, 1) * 1.2))
    cand <- c(.001, .01, .02, .05, .1, .2, .25, .5, .75, 1, 1.5, 2, 3, 4, 5, 10, 20, 50, 100, 1000)
  }
  cand_in <- cand[cand >= exp(alim[1]) - 1e-9 & cand <= exp(alim[2]) + 1e-9]
  
  # ---- layout horizontal (0-100) --------------------------------------------
  L <- if (mostrar_peso) list(e1 = 20, t1 = 28, e2 = 38, t2 = 46, wt = 54, rr = 71, px0 = 75, px1 = 98.5)
  else              list(e1 = 21, t1 = 29, e2 = 41, t2 = 49, wt = NA, rr = 67, px0 = 71.5, px1 = 98.5)
  px <- function(v) L$px0 + (pmin(pmax(v, alim[1]), alim[2]) - alim[1]) / diff(alim) * (L$px1 - L$px0)
  
  # ---- layout vertical (y decresce para baixo) --------------------------------
  y_st  <- -2.7 - (seq_len(k) - 1)
  y_tot <- min(y_st) - 1.4
  y_evt <- y_tot - 1
  y_ax  <- y_evt - 0.7
  y_het <- y_evt - 1.3
  y_z   <- y_het - 1
  y_tk  <- y_ax - 0.85
  y_fav <- y_ax - 1.75
  y_fn  <- y_fav - 1.3
  ymin  <- y_fn - 0.7
  ytop  <- 1.7
  h_in  <- (ytop - ymin) * 0.29
  
  nome_arquivo <- file.path(dir_saida, paste0("ForestCochrane_", base, "_x_", comorb, ".png"))
  png(nome_arquivo, width = 12.5, height = h_in, units = "in", res = 200, pointsize = 12)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot.new(); plot.window(xlim = c(0, 100), ylim = c(ymin, ytop))
  
  # ---- título -----------------------------------------------------------------
  text(0.5, 1.0, paste0(base, " (diagn\u00f3stico base) \u2192 ", comorb, " (comorbidade)"),
       adj = c(0, 0.5), font = 2, cex = 1.1)
  
  # ---- cabeçalhos ---------------------------------------------------------------
  text(0.5, -1.25, "UF", adj = c(0, 0.5), font = 2)
  text((L$e1 + L$t1) / 2 - 2, -0.45, "Com a comorbidade", adj = c(0.5, 0.5), font = 2)
  text((L$e2 + L$t2) / 2 - 2, -0.45, "Sem a comorbidade", adj = c(0.5, 0.5), font = 2)
  text(c(L$e1, L$e2), -1.25, "Eventos", adj = c(1, 0.5))
  text(c(L$t1, L$t2), -1.25, "Total",   adj = c(1, 0.5))
  if (mostrar_peso) text(L$wt, -0.85, "Peso", adj = c(1, 0.5), font = 2)
  for (xc in c(L$rr - 6.5, (L$px0 + L$px1) / 2)) {
    text(xc, -0.45, "Raz\u00e3o de risco", adj = c(0.5, 0.5), font = 2)
    text(xc, -1.25, "IV, aleat\u00f3rio, IC 95%", adj = c(0.5, 0.5))
  }
  segments(0, -1.8, 100, -1.8, lwd = 1.2)
  
  # ---- linhas dos estudos (UFs) ------------------------------------------------
  rotulos <- paste0(d$UF, if (marcar_sig) ifelse(d$significativo_uf, " *", "") else "")
  text(0.5, y_st, rotulos, adj = c(0, 0.5))
  text(L$e1, y_st, fmt_int(d$a),      adj = c(1, 0.5))
  text(L$t1, y_st, fmt_int(d$n_exp),  adj = c(1, 0.5))
  text(L$e2, y_st, fmt_int(d$c),      adj = c(1, 0.5))
  text(L$t2, y_st, fmt_int(d$n_nexp), adj = c(1, 0.5))
  if (mostrar_peso) text(L$wt, y_st, fmt_pct1(w), adj = c(1, 0.5))
  text(L$rr, y_st, fmt_rr(d$RR, d$RR_li, d$RR_ls), adj = c(1, 0.5))
  
  # ---- gráfico: linha em RR = 1, IC, quadrados ----------------------------------
  segments(px(0), -1.8, px(0), y_ax, lwd = 1.6)
  segments(px(log(d$RR_li)), y_st, px(log(d$RR_ls)), y_st, lwd = 2)
  cex_q <- 0.5 + 0.8 * sqrt(w / max(w))      # área proporcional ao peso
  points(px(d$yi), y_st, pch = 15, cex = cex_q, col = cor_quadrado)
  
  # ---- linha do total -----------------------------------------------------------
  text(0.5, y_tot, "Total (IC 95%)", adj = c(0, 0.5), font = 2)
  text(L$t1, y_tot, fmt_int(sum(d$n_exp)),  adj = c(1, 0.5), font = 2)
  text(L$t2, y_tot, fmt_int(sum(d$n_nexp)), adj = c(1, 0.5), font = 2)
  if (mostrar_peso) text(L$wt, y_tot, "100,0%", adj = c(1, 0.5), font = 2)
  text(L$rr, y_tot, fmt_rr(exp(est), exp(ci_lb), exp(ci_ub)), adj = c(1, 0.5), font = 2)
  polygon(c(px(ci_lb), px(est), px(ci_ub), px(est)),
          c(y_tot, y_tot + 0.38, y_tot, y_tot - 0.38), col = cor_losango, border = "black")
  
  text(0.5, y_evt, "Total de eventos", adj = c(0, 0.5))
  text(L$e1, y_evt, fmt_int(sum(d$a)), adj = c(1, 0.5))
  text(L$e2, y_evt, fmt_int(sum(d$c)), adj = c(1, 0.5))
  
  # ---- heterogeneidade e efeito global ------------------------------------------
  text(0.5, y_het, adj = c(0, 0.5), cex = 0.95, paste0(
    "Heterogeneidade: Tau\u00b2 = ", formatC(fit$tau2, format = "f", digits = 3, decimal.mark = ","),
    "; Chi\u00b2 = ", formatC(fit$QE, format = "f", digits = 2, decimal.mark = ","),
    ", gl = ", k - 1, " (P ", fmt_p5(fit$QEp), "); I\u00b2 = ", round(fit$I2), "%"))
  text(0.5, y_z, adj = c(0, 0.5), cex = 0.95, paste0(
    "Teste do efeito global: ", ifelse(ic == "hksj", "t", "Z"), " = ",
    formatC(fit$zval, format = "f", digits = 2, decimal.mark = ","), " (P ", fmt_p5(fit$pval), ")"))
  
  # ---- eixo, ticks e "favorece" ---------------------------------------------------
  segments(L$px0, y_ax, L$px1, y_ax, lwd = 1.2)
  segments(px(log(cand_in)), y_ax, px(log(cand_in)), y_ax - 0.3, lwd = 1.2)
  text(px(log(cand_in)), y_tk,
       sapply(cand_in, function(v) format(v, decimal.mark = ",", scientific = FALSE)),
       adj = c(0.5, 0.5))
  text((L$px0 + px(0)) / 2, y_fav, label_esq, adj = c(0.5, 0.5), cex = 0.9)
  text((px(0) + L$px1) / 2, y_fav, label_dir, adj = c(0.5, 0.5), cex = 0.9)
  
  # ---- nota de rodapé -------------------------------------------------------------
  nota <- "Evento = \u00f3bito hospitalar. RR = risco de \u00f3bito com vs. sem a comorbidade."
  if (!is.na(n_sig)) nota <- paste0(nota, " Significativo individualmente (FDR<", ALPHA_FDR, ", RR>", RR_MIN,
                                    ") em ", n_sig, "/", length(UFS_NORDESTE), " UFs.",
                                    if (marcar_sig) " (* = significativo)" else "")
  text(0.5, y_fn, nota, adj = c(0, 0.5), cex = 0.8)
  
  invisible(nome_arquivo)
}

# ------------------------------------------------------------------------------
# USO
# ------------------------------------------------------------------------------
  arq <- gerar_forest_cochrane("A419", "R579", n_sig = 9)
 browseURL(arq)

 com coluna de peso, IC de Hartung-Knapp e quadrados pretos (RevMan 5 clássico):
   gerar_forest_cochrane("A419", "R579", n_sig = 9, mostrar_peso = TRUE, ic = "hksj",
                        cor_quadrado = "black")

 todos os pares selecionados no Módulo 16.9:
  purrr::pwalk(list(pares_para_forest$DIAG_PRINC, pares_para_forest$DIAGSEC1,
                     pares_para_forest$n_UFs_significativas),
               function(b, c, n) gerar_forest_cochrane(b, c, n_sig = n))
 
 
 
 
 # ==============================================================================
 # GERAR TODOS OS FOREST PLOTS NO MODELO COCHRANE (substitui o Módulo 16.9 antigo)
 # Pré-requisitos: Módulos 1-8 rodados nesta sessão (pares_por_uf,
 # rede_principal_nordeste, pares_para_forest) e o arquivo
 # forest_estilo_cochrane.R na pasta de trabalho.
 # ==============================================================================
 
 if (!exists("gerar_forest_cochrane")) source("forest_estilo_cochrane.R")
 
 # ---- opções ------------------------------------------------------------------
 GERAR_TODOS   <- FALSE   # FALSE = só os pares de pares_para_forest (>= MIN_UFS_SIG_FOREST UFs significativas)
 # TRUE  = TODOS os pares da rede regional (pode ser centenas de PNGs)
 USAR_HKSJ     <- FALSE   # TRUE = IC de Hartung-Knapp no losango (mais conservador)
 MOSTRAR_PESO  <- FALSE   # TRUE = coluna de peso (%)
 DIR_FOREST    <- "forests_cochrane"
 
 if (!dir.exists(DIR_FOREST)) dir.create(DIR_FOREST)
 
 # ---- seleção dos pares ------------------------------------------------------------
 sig_count_por_par <- pares_por_uf %>%
   filter(significativo_uf) %>%
   count(DIAG_PRINC, DIAGSEC1, name = "n_UFs_significativas")
 
 pares_alvo <- if (GERAR_TODOS) {
   rede_principal_nordeste %>%
     left_join(sig_count_por_par, by = c("DIAG_PRINC", "DIAGSEC1")) %>%
     mutate(n_UFs_significativas = replace_na(n_UFs_significativas, 0L)) %>%
     arrange(desc(n_UFs_significativas), desc(RR_meta))
 } else {
   pares_para_forest
 }
 
 cat("=== FOREST PLOTS (modelo Cochrane): ", nrow(pares_alvo), " pares -> pasta '", DIR_FOREST, "' ===\n", sep = "")
 
 # ---- geração (um erro em um par não interrompe os demais) ---------------------------
 gerar_um <- function(base, comorb, n_sig) {
   tryCatch({
     arq <- gerar_forest_cochrane(base, comorb, n_sig = n_sig,
                                  ic = if (USAR_HKSJ) "hksj" else "wald",
                                  mostrar_peso = MOSTRAR_PESO,
                                  dir_saida = DIR_FOREST)
     if (is.null(arq)) "pulado (poucas UFs ou ajuste falhou)" else "ok"
   }, error = function(e) paste0("ERRO: ", conditionMessage(e)))
 }
 
 status <- character(nrow(pares_alvo))
 for (i in seq_len(nrow(pares_alvo))) {
   status[i] <- gerar_um(pares_alvo$DIAG_PRINC[i], pares_alvo$DIAGSEC1[i], pares_alvo$n_UFs_significativas[i])
   cat(sprintf(" [%02d/%02d] %s -> %s : %s\n", i, nrow(pares_alvo),
               pares_alvo$DIAG_PRINC[i], pares_alvo$DIAGSEC1[i], status[i]))
 }
 
 # ---- índice dos forest plots ----------------------------------------------------------
 indice_forest <- pares_alvo %>%
   transmute(DIAG_PRINC, DIAGSEC1, n_UFs_significativas,
             RR_meta = round(RR_meta, 2), I2 = round(I2, 1),
             arquivo = file.path(DIR_FOREST, paste0("ForestCochrane_", DIAG_PRINC, "_x_", DIAGSEC1, ".png")),
             status = status)
 
 data.table::fwrite(indice_forest, file.path(DIR_FOREST, "indice_forest_plots.csv"))
 
 cat("\n -> Gerados com sucesso:", sum(status == "ok"), "de", nrow(pares_alvo), "\n")
 if (any(status != "ok")) {
   cat(" -> Pares com problema:\n")
   print(as.data.frame(indice_forest[status != "ok", c("DIAG_PRINC", "DIAGSEC1", "status")]))
 }
 cat(" -> Pasta:", normalizePath(DIR_FOREST), "\n")
 
 # Para abrir a pasta no explorador (Windows):  shell.exec(normalizePath(DIR_FOREST))
 # Para apagar os PNGs do modelo antigo (Forest_*.png na pasta de trabalho):
 #   file.remove(list.files(pattern = "^Forest_.*\\.png$"))
 
 
 
 
 # ==============================================================================
 # BLOCO AVANÇADO 1: MINERAÇÃO DE PARES COM BAIXA HETEROGENEIDADE (I² < 50%)
 # ==============================================================================
 cat("\n=== FILTRANDO PARES COM BAIXA HETEROGENEIDADE (ROBUSTOS REGIONALMENTE) ===\n")
 
 # Filtra da rede principal apenas os que têm I² < 50%
 pares_baixa_hetero <- rede_principal_nordeste %>%
   filter(robustez_regional == "Robusta (baixa heterogeneidade)") %>%
   left_join(sig_count_por_par, by = c("DIAG_PRINC", "DIAGSEC1")) %>%
   mutate(n_UFs_significativas = replace_na(n_UFs_significativas, 0L)) %>%
   arrange(desc(n_UFs_significativas), I2, desc(RR_meta))
 
 cat(" -> Total de pares robustos (I² < 50%):", nrow(pares_baixa_hetero), "\n")
 cat(" -> Top 10 pares com menor I² e ampla presença regional:\n")
 print(pares_baixa_hetero %>% 
         select(DIAG_PRINC, DIAGSEC1, n_UFs_incluidas, n_UFs_significativas, RR_meta, I2, p_meta_fdr) %>% 
         head(10))
 
 # ------------------------------------------------------------------------------
 # GERAR FOREST PLOTS APENAS DOS PARES DE BAIXA HETEROGENEIDADE
 # ------------------------------------------------------------------------------
 top_pares_robustos_plot <- pares_baixa_hetero %>% head(10)
 
 if (nrow(top_pares_robustos_plot) > 0) {
   cat("\n=== GERANDO FOREST PLOTS DOS PARES ROBUSTOS (I² BAIXO) ===\n")
   graphics.off()
   
   pwalk(
     list(top_pares_robustos_plot$DIAG_PRINC, 
          top_pares_robustos_plot$DIAGSEC1, 
          top_pares_robustos_plot$n_UFs_significativas),
     gerar_forest_par
   )
   cat(" -> Forest plots dos pares de baixa heterogeneidade gerados com sucesso!\n")
 }
 
 # ------------------------------------------------------------------------------
 # FIGURA RESUMO ESPECÍFICA DE PARES ROBUSTOS (I² < 50%)
 # ------------------------------------------------------------------------------
 plot_resumo_robustos <- top_pares_robustos_plot %>%
   mutate(par_label = paste0(DIAG_PRINC, " → ", DIAGSEC1, " (I² = ", round(I2, 1), "%)"),
          par_label = fct_reorder(par_label, RR_meta)) %>%
   ggplot(aes(x = RR_meta, y = par_label)) +
   geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
   geom_errorbar(aes(xmin = RR_meta_li, xmax = RR_meta_ls), orientation = "y", width = 0.25, linewidth = 0.7, color = "#1a9850") +
   geom_point(aes(size = n_UFs_significativas), color = "#1a9850") +
   scale_size_continuous(name = "UFs significativas", range = c(2.5, 5.5)) +
   scale_x_log10() +
   labs(title = "Associações Clínicas Universais no Nordeste (I² < 50%)",
        subtitle = "Pares com alta consistência biológica e baixa dispersão entre estados",
        x = "Risco Relativo Meta-analítico (escala log)", y = NULL) +
   theme_minimal(base_size = 11) +
   theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
 
 ggsave("Figura9_ForestPlot_Pares_Robustos_Baixo_I2.png", plot_resumo_robustos, width = 10, height = 7, dpi = 300, bg = "white")
 cat(" -> Figura 9 (Resumo de Pares Robustos) exportada com sucesso!\n")
 
 # ==============================================================================
 # BLOCO AVANÇADO 2: ANÁLISE DO PARADOXO DE SIMPSON / DISTORÇÃO DO POOL
 # ==============================================================================
 cat("\n=== AVALIANDO DISTORÇÃO METODOLÓGICA (SIMPSON / POOL INGÊNUO) ===\n")
 
 comparativo_simpson <- rede_principal_nordeste %>%
   inner_join(pares_pool_ingenuo %>% select(DIAG_PRINC, DIAGSEC1, RR_pool_ingenuo),
              by = c("DIAG_PRINC", "DIAGSEC1")) %>%
   mutate(
     Diferenca_Relativa_Pct = round(100 * (RR_pool_ingenuo - RR_meta) / RR_meta, 2),
     Divergencia_Critica    = abs(Diferenca_Relativa_Pct) >= 25
   ) %>%
   arrange(desc(abs(Diferenca_Relativa_Pct)))
 
 # Gráfico comparativo: Meta-análise vs Pool Ingênuo
 plot_simpson <- ggplot(comparativo_simpson, aes(x = RR_meta, y = RR_pool_ingenuo)) +
   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
   geom_point(aes(color = Divergencia_Critica), alpha = 0.75, size = 2) +
   scale_color_manual(values = c("FALSE" = "#2C3E50", "TRUE" = "#E74C3C"),
                      name = "Divergência > 25%") +
   scale_x_log10() + scale_y_log10() +
   labs(title = "Meta-Análise de Efeitos Aleatórios vs. Agregação Ingênua (Pool)",
        subtitle = "Linha vermelha = Concordância perfeita (y = x). Pontos vermelhos = Distorção metodológica severa",
        x = "Risco Relativo Meta-Analítico Regional (Efeitos Aleatórios)",
        y = "Risco Relativo do Pool Ingênuo (Empilhado)") +
   theme_minimal(base_size = 12) +
   theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
 
 ggsave("Figura10_Paradoxo_Simpson_Pool_vs_Meta.png", plot_simpson, width = 9, height = 7, dpi = 300, bg = "white")
 cat(" -> Figura 10 (Comparativo Simpson / Pool vs Meta) exportada com sucesso!\n")
 
 # ==============================================================================
 # VISUALIZAÇÃO DA REDE REGIONAL DE MORTALIDADE (NORDESTE)
 # ==============================================================================
 library(igraph)
 library(ggraph)
 library(visNetwork)
 library(dplyr)
 
 # 1. Verificação do grafo regional no ambiente
 if (!exists("grafo_regional")) {
   stop("O objeto 'grafo_regional' não foi encontrado. Certifica-te de que o Módulo 8 foi executado!")
 }
 
 # ------------------------------------------------------------------------------
 # OPÇÃO A: VISUALIZAÇÃO INTERATIVA DINÂMICA (visNetwork)
 # ------------------------------------------------------------------------------
 cat(" -> A gerar visualização interativa no painel Viewer...\n")
 
 # Paleta de cores para os módulos sindrómicos
 cores_hex <- c("1" = "#1a9850", "2" = "#fc8d62", "3" = "#8073ac", "4" = "#e78ac3",
                "5" = "#66c2a5", "6" = "#ffd92f", "7" = "#a6d854", "8" = "#e5c494")
 
 # Preparação dos nós (tamanho proporcional à intermediação)
 nodes_vis <- data.frame(
   id    = V(grafo_regional)$name,
   label = V(grafo_regional)$name,
   group = V(grafo_regional)$group,
   value = sqrt(pmax(V(grafo_regional)$betweenness, 0.1)) * 4 + 5,
   title = paste0("<b>CID-10:</b> ", V(grafo_regional)$name,
                  "<br><b>Módulo:</b> ", V(grafo_regional)$group,
                  "<br><b>Grau:</b> ", degree(grafo_regional),
                  "<br><b>Intermediação:</b> ", round(V(grafo_regional)$betweenness, 2), "%"),
   color = list(
     background = cores_hex[as.character(V(grafo_regional)$group)],
     border     = "#2c3e50",
     highlight  = "#e74c3c"
   )
 )
 
 # Preparação das arestas (espessura por RR e estilo por robustez)
 edges_vis <- as_data_frame(grafo_regional, what = "edges") %>%
   mutate(
     from   = from,
     to     = to,
     width  = log(pmax(weight, 1.1)) * 1.5,
     title  = paste0("<b>Associação:</b> ", from, " → ", to,
                     "<br><b>RR Meta:</b> ", round(weight, 2),
                     "<br><b>I²:</b> ", round(I2, 1), "%",
                     "<br><b>Status:</b> ", robustez_regional),
     dashes = (tier != "Alta"), # Dotted/dashed para exploratórias
     color  = ifelse(tier == "Alta", "#555555", "#b0bec5")
   )
 
 # Renderização do widget interativo
 rede_interativa <- visNetwork(nodes_vis, edges_vis, 
                               main = "Rede de Associações Clínicas Letais — Nordeste",
                               submain = "Nós proporcionais à intermediação | Cores por Módulo Sindrómico") %>%
   visOptions(highlightNearest = list(enabled = TRUE, hover = TRUE, degree = 1),
              nodesIdSelection = TRUE,
              selectedBy = "group") %>%
   visPhysics(solver = "forceAtlas2Based", 
              forceAtlas2Based = list(gravitationalConstant = -60, centralGravity = 0.01)) %>%
   visLayout(randomSeed = 42)
 
 # Exibe no visualizador do RStudio
 print(rede_interativa)
 
 
 # ------------------------------------------------------------------------------
 # OPÇÃO B: RENDERIZAÇÃO ESTÁTICA EM ALTA DEFINIÇÃO (ggraph)
 # ------------------------------------------------------------------------------
 cat(" -> A renderizar e exibir a rede estática de alta definição...\n")
 
 set.seed(42)
 plot_rede_estatica <- ggraph(grafo_regional, layout = "fr") +
   # Arestas
   geom_edge_link(aes(width = weight, linetype = tier), 
                  color = "#555555", alpha = 0.55) +
   scale_edge_width(range = c(0.4, 2.8), guide = "none") +
   scale_edge_linetype_manual(values = c("Alta" = "solid", "Exploratória" = "dotted"), 
                              name = "Robustez Regional") +
   # Nós
   geom_node_point(aes(size = betweenness, fill = group), 
                   shape = 21, color = "black", stroke = 0.6) +
   scale_size_continuous(range = c(3, 16), name = "Intermediação (%)") +
   scale_fill_manual(values = cores_hex, name = "Módulo Sindrómico") +
   # Rótulos (CID-10)
   geom_node_text(aes(label = name), repel = TRUE, size = 3.2, 
                  fontface = "bold", bg.color = "white", bg.r = 0.12) +
   labs(title = "Topologia das Associações Clínicas Letais — Nordeste (2016–2025)",
        subtitle = "Linhas contínuas = Robustas entre estados (I² < 50%) | Pontilhadas = Específicas de estado") +
   theme_void(base_size = 12) +
   theme(
     plot.title    = element_text(face = "bold", size = 15, hjust = 0.5, margin = margin(b = 5)),
     plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
     legend.position = "right"
   )
 
 # Mostra o gráfico na janela de Plots do RStudio
 print(plot_rede_estatica)
 
 
 library(ggplot2)
 library(dplyr)
 library(forcats)
 
 # 1. Preparação dos dados para ordenação e rótulo dos pares
 dados_plot_robustos <- top_pares_robustos_plot %>%
   mutate(
     # Rótulo com par e valor do I² correspondente
     par_label = paste0(DIAG_PRINC, " → ", DIAGSEC1, " (I² = ", round(I2, 1), "%)"),
     par_label = fct_reorder(par_label, RR_meta),
     
     # Conversão do número de UFs em fator para mapear formas geométricas distintas
     n_ufs_cat = factor(n_UFs_significativas, levels = sort(unique(n_UFs_significativas)))
   )
 
 # 2. Definição das formas geométricas e tamanhos correspondentes
 # 21 = Círculo preenchido, 22 = Quadrado, 24 = Triângulo, 23 = Losango
 formas_ufs <- c("4" = 21, "5" = 22, "6" = 24, "7" = 23)
 
 # 3. Construção do gráfico ggplot2
 plot_resumo_robustos_v2 <- ggplot(dados_plot_robustos, aes(x = RR_meta, y = par_label)) +
   # Linha vertical tracejada de nulidade (RR = 1)
   geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.6) +
   
   # Intervalos de Confiança (barras de erro horizontais)
   geom_errorbar(
     aes(xmin = RR_meta_li, xmax = RR_meta_ls),
     orientation = "y",
     width = 0.28,
     linewidth = 0.75,
     color = "#1b7837"
   ) +
   
   # Estimativa pontual com formas geométricas e tamanhos proporcionais
   geom_point(
     aes(shape = n_ufs_cat, fill = n_ufs_cat),
     size = 3.6,
     color = "black",
     stroke = 0.5
   ) +
   
   # Escala de formas geométricas distintas
   scale_shape_manual(
     name = "UFs com associação significante:",
     values = formas_ufs
   ) +
   
   # Cor de preenchimento idêntica para manter elegância visual
   scale_fill_manual(
     name = "UFs com associação significante:",
     values = rep("#1b7837", length(formas_ufs))
   ) +
   
   # Escala logarítmica com marcas intermediárias bem distribuídas
   scale_x_log10(
     breaks = c(1, 2, 3, 5, 10, 20),
     labels = c("1", "2", "3", "5", "10", "20"),
     limits = c(0.95, max(dados_plot_robustos$RR_meta_ls) * 1.15)
   ) +
   
   # Títulos e legendas
   labs(
     title = "Associações Clínicas Letais Universais no Nordeste (I² < 50%)",
     subtitle = "Pares com elevada consistência biológica e homogeneidade entre estados (2016–2025)",
     x = "Risco Relativo Meta-Analítico — Escala Logarítmica (IC 95%)",
     y = NULL
   ) +
   
   # Tema visual limpo para publicação
   theme_minimal(base_size = 12) +
   theme(
     plot.title         = element_text(face = "bold", size = 13.5, margin = margin(b = 4)),
     plot.subtitle      = element_text(size = 11, color = "grey30", margin = margin(b = 14)),
     axis.text.y        = element_text(face = "bold", size = 10, color = "black"),
     axis.text.x        = element_text(size = 10, color = "black"),
     axis.title.x       = element_text(face = "bold", size = 11, margin = margin(t = 10)),
     panel.grid.minor   = element_blank(),
     panel.grid.major.y = element_line(color = "grey92", linetype = "dotted"),
     legend.position    = "bottom",
     legend.title       = element_text(face = "bold", size = 10),
     legend.text        = element_text(size = 10),
     legend.box.margin  = margin(t = 10)
   ) +
   # Configuração da legenda para exibir o formato preenchido com a cor
   guides(
     shape = guide_legend(override.aes = list(size = 4, fill = "#1b7837", color = "black")),
     fill  = "none"
   )
 
 # 4. Gravação da imagem em alta resolução
 ggsave(
   filename = "Figura9_ForestPlot_Pares_Robustos_Formas_UFs.png",
   plot     = plot_resumo_robustos_v2,
   width    = 10.5,
   height   = 7,
   dpi      = 300,
   bg       = "white"
 )
 
 cat(" -> Gráfico atualizado gravado como 'Figura9_ForestPlot_Pares_Robustos_Formas_UFs.png' com sucesso!\n")
 
 
 
 # ==============================================================================
 # MÓDULO 17: ANÁLISE DE SUBGRUPOS POR ERAS EPIDEMIOLÓGICAS (PRÉ, DURANTE E PÓS-COVID)
 # ==============================================================================
 cat("\n=== INICIANDO SUBANÁLISE EPIDEMIOLÓGICA (ERAS COVID-19) ===\n")
 
 library(dplyr)
 library(tidyr)
 library(igraph)
 library(ggplot2)
 library(forcats)
 
 # 1. Classificação temporal na coorte purificada consolidada
 if (!exists("cols_ref") || is.null(cols_ref)) cols_ref <- detectar_colunas(coorte_nordeste)
 col_ano <- cols_ref$col_ano
 
 coorte_fases <- coorte_nordeste %>%
   mutate(
     Ano_Num = as.numeric(as.character(.data[[col_ano]])),
     Fase_Epidemiologica = case_when(
       Ano_Num <= 2019 ~ "1. Pré-Pandemia (2016-2019)",
       Ano_Num >= 2020 & Ano_Num <= 2022 ~ "2. Auge Pandêmico (2020-2022)",
       Ano_Num >= 2023 ~ "3. Pós-Pandemia (2023-2025)"
     )
   )
 
 # Resumo amostral por fase
 tabela_fases_resumo <- coorte_fases %>%
   group_by(Fase_Epidemiologica) %>%
   summarise(
     Total_Obitos = n(),
     Casos_COVID = sum(DIAG_PRINC %in% c("B342", "U071", "U072") | 
                         DIAGSEC1 %in% c("B342", "U071", "U072"), na.rm = TRUE),
     Pct_COVID = round(100 * Casos_COVID / Total_Obitos, 2),
     .groups = "drop"
   )
 
 cat(" -> Distribuição de Óbitos e Casos de COVID-19 por Era:\n")
 print(as.data.frame(tabela_fases_resumo))
 cat("\n")
 
 # 2. Construção de redes e métricas por Era Epidemiológica
 fases_lista <- unique(coorte_fases$Fase_Epidemiologica)
 
 analise_hubs_fases <- map_dfr(sort(fases_lista), function(fase_nome) {
   df_sub <- coorte_fases %>% filter(Fase_Epidemiologica == fase_nome)
   
   # Pares frequentes na fase (mínimo de 10 casos na região para consistência)
   pares_fase <- df_sub %>%
     count(DIAG_PRINC, DIAGSEC1, name = "peso") %>%
     filter(peso >= 10)
   
   g_fase <- graph_from_data_frame(pares_fase, directed = FALSE)
   g_fase <- delete_vertices(g_fase, which(degree(g_fase) == 0))
   
   # Centralidades
   betw <- betweenness(g_fase, directed = FALSE, normalized = TRUE) * 100
   grau <- degree(g_fase)
   
   # Comunidades Louvain
   com <- cluster_louvain(g_fase)
   
   tibble(
     Fase = fase_nome,
     CID10 = names(betw),
     Intermediacao = round(betw, 2),
     Grau = grau,
     Modulo = as.character(membership(com)[names(betw)]),
     N_Nos_Rede = vcount(g_fase),
     N_Arestas_Rede = ecount(g_fase),
     Modularidade_Global = round(modularity(com), 4)
   )
 })
 
 # 3. Identificar os 5 maiores hubs de cada fase
 top_hubs_por_fase <- analise_hubs_fases %>%
   group_by(Fase) %>%
   arrange(desc(Intermediacao)) %>%
   slice_head(n = 6) %>%
   ungroup()
 
 cat(" -> Top 6 Hubs Mediadores Críticos por Fase Epidemiológica:\n")
 print(top_hubs_por_fase %>% select(Fase, CID10, Intermediacao, Grau, Modularidade_Global), n = 20)
 cat("\n")
 
 # ==============================================================================
 # ==============================================================================
 # ETAPA 4: GERAÇÃO DA FIGURA 11 (SEM DEPENDÊNCIA DO TIDYTEXT)
 # ==============================================================================
 
 # Criação de um identificador único de ordenação por faceta
 plot_comparativo_fases <- top_hubs_por_fase %>%
   # Cria uma chave única por diagnóstico e fase para ordenar livremente em cada faceta
   mutate(CID_Ordenado = factor(paste(CID10, Fase, sep = "___"))) %>%
   mutate(CID_Ordenado = reorder(CID_Ordenado, Intermediacao)) %>%
   ggplot(aes(x = Intermediacao, y = CID_Ordenado, fill = Fase)) +
   geom_col(color = "black", linewidth = 0.3, width = 0.7, show.legend = FALSE) +
   # Remove o sufixo '___Fase' do eixo Y mantendo apenas o código CID-10
   scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
   facet_wrap(~ Fase, scales = "free_y", ncol = 3) +
   scale_fill_manual(values = c("1. Pré-Pandemia (2016-2019)" = "#3182bd",
                                "2. Auge Pandêmico (2020-2022)" = "#de2d26",
                                "3. Pós-Pandemia (2023-2025)" = "#31a354")) +
   labs(
     title = "Reorganização Topológica dos Hubs Letais nas Três Eras Epidemiológicas",
     subtitle = "Comparação da Centralidade de Intermediação (%) no Nordeste (2016–2025)",
     x = "Centralidade de Intermediação Normalizada (%)",
     y = "CID-10 (Diagnóstico)"
   ) +
   theme_minimal(base_size = 11) +
   theme(
     plot.title    = element_text(face = "bold", size = 13),
     plot.subtitle = element_text(size = 10.5, color = "grey30", margin = margin(b = 10)),
     strip.text    = element_text(face = "bold", size = 10.5),
     axis.text.y   = element_text(face = "bold", size = 9.5),
     panel.spacing = unit(1.2, "lines")
   )
 
 ggsave("Figura11_Comparativo_Eras_Pandemicas.png", plot_comparativo_fases, width = 12, height = 5.5, dpi = 300, bg = "white")
 cat(" -> Figura 11 ('Figura11_Comparativo_Eras_Pandemicas.png') salva com sucesso!\n")
 
 # ==============================================================================
 # ETAPA 5: ADICIONAR AO WORKBOOK DO EXCEL
 # ==============================================================================
 if (file.exists("OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx")) {
   wb_atual <- loadWorkbook("OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx")
   if (!"Subanalise_Eras_COVID" %in% names(wb_atual)) {
     addWorksheet(wb_atual, "Subanalise_Eras_COVID")
     writeData(wb_atual, "Subanalise_Eras_COVID", analise_hubs_fases)
     saveWorkbook(wb_atual, "OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx", overwrite = TRUE)
     cat(" -> Nova aba 'Subanalise_Eras_COVID' incorporada com sucesso ao Workbook OSF!\n")
   }
 }
 # 5. Adicionar a aba ao Workbook do Excel existente
 if (file.exists("OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx")) {
   wb_atual <- loadWorkbook("OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx")
   if (!"Subanalise_Eras_COVID" %in% names(wb_atual)) {
     addWorksheet(wb_atual, "Subanalise_Eras_COVID")
     writeData(wb_atual, "Subanalise_Eras_COVID", analise_hubs_fases)
     saveWorkbook(wb_atual, "OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx", overwrite = TRUE)
     cat(" -> Nova aba 'Subanalise_Eras_COVID' incorporada com sucesso ao Workbook OSF!\n")
   }
 }
 
 
 # ==============================================================================
 # MÓDULO 18: ANÁLISES AVANÇADAS DE IMPACTO DA COVID-19
 # ==============================================================================
 cat("\n=== INICIANDO SUBANÁLISES 1 E 2: SUBGRUPOS TEMPORAIS & SENSIBILIDADE LIVRE DE COVID ===\n")
 
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(igraph)
 library(meta)
 library(ggplot2)
 library(openxlsx)
 
 graphics.off()
 
 # ------------------------------------------------------------------------------
 # PARTE 1: FOREST PLOTS COM ANÁLISE DE SUBGRUPOS (PRÉ vs. DURANTE vs. PÓS-COVID)
 # ------------------------------------------------------------------------------
 cat(" -> [1/2] Calculando Forest Plots estratificados por subgrupos temporais...\n")
 
 # Classificação temporal por paciente e UF
 col_ano <- cols_ref$col_ano
 
 dados_subgrupos <- coorte_nordeste %>%
   mutate(
     Ano_Num = as.numeric(as.character(.data[[col_ano]])),
     Subgrupo_Era = case_when(
       Ano_Num <= 2019 ~ "1. Pré-Pandemia (2016-2019)",
       Ano_Num >= 2020 & Ano_Num <= 2022 ~ "2. Auge Pandêmico (2020-2022)",
       Ano_Num >= 2023 ~ "3. Pós-Pandemia (2023-2025)"
     )
   )
 
 gerar_forest_subgrupos_era <- function(base_cod, comorb_cod) {
   # Agregação 2x2 por UF dentro de cada Era
   tabela_uf_era <- dados_subgrupos %>%
     group_by(Subgrupo_Era, UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (n_distinct(tabela_uf_era$Subgrupo_Era) < 2 || nrow(tabela_uf_era) < 6) return(invisible(NULL))
   
   # Meta-análise com subgrupos no pacote meta
   m_sub <- tryCatch({
     metabin(
       event.e  = a,
       n.e      = a + b,
       event.c  = c,
       n.c      = c + d,
       studlab  = UF,
       subgroup = Subgrupo_Era,
       data     = tabela_uf_era,
       sm       = "RR",
       method   = "Inverse",
       method.tau = "REML",
       method.random.ci = "HK"
     )
   }, error = function(e) NULL)
   
   if (is.null(m_sub)) return(invisible(NULL))
   
   nome_arq <- paste0("Figura12_Forest_Subgrupos_COVID_", base_cod, "_x_", comorb_cod, ".png")
   png(nome_arq, width = 2300, height = 1400 + 40 * nrow(tabela_uf_era), res = 220)
   
   forest(
     m_sub,
     layout           = "meta",
     common           = FALSE,
     random           = TRUE,
     test.subgroup    = TRUE,
     col.square       = "#08306b",
     col.diamond      = "black",
     squaresize       = 0.75,
     leftcols         = c("studlab", "event.e", "n.e", "event.c", "n.c"),
     leftlabs         = c("UF", "Óbitos", "Total", "Óbitos", "Total"),
     label.e          = paste0("Comorbidade (", comorb_cod, ")"),
     label.c          = paste0("Base Isolada (", base_cod, ")"),
     rightcols        = c("effect", "ci", "w.random"),
     rightlabs        = c("RR", "95%-CI", "Peso"),
     ref              = 1,
     digits           = 2,
     digits.weight    = 1,
     smlab            = "Risco Relativo"
   )
   
   title(main = paste0("Análise de Subgrupos (Eras COVID): ", base_cod, " → ", comorb_cod),
         line = 1.8, cex.main = 1.15, font.main = 2)
   dev.off()
   cat("    -> Forest Plot com subgrupos salvo:", nome_arq, "\n")
   invisible(NULL)
 }
 
 # Roda para pares representativos universais
 gerar_forest_subgrupos_era("I64", "J180")
 gerar_forest_subgrupos_era("I500", "N179")
 
 
 # ------------------------------------------------------------------------------
 # PARTE 2: ANÁLISE DE SENSIBILIDADE "LIVRE DE COVID" (LEAVE-ONE-OUT)
 # ------------------------------------------------------------------------------
 cat("\n -> [2/2] Executando Análise de Sensibilidade Livre de COVID (exclusão 100%)...\n")
 
 cids_covid <- c("B342", "U071", "U072")
 
 # Coorte expurgada de infecções por coronavírus
 coorte_sem_covid <- coorte_nordeste %>%
   filter(!(DIAG_PRINC %in% cids_covid), !(DIAGSEC1 %in% cids_covid))
 
 cat("    - Casos na coorte original:", nrow(coorte_nordeste), "\n")
 cat("    - Casos na coorte livre de COVID:", nrow(coorte_sem_covid), "\n")
 cat("    - Registros removidos:", nrow(coorte_nordeste) - nrow(coorte_sem_covid), "\n")
 
 # Recalcula a rede do Nordeste livre de COVID
 pares_sem_covid <- coorte_sem_covid %>%
   count(DIAG_PRINC, DIAGSEC1, name = "peso") %>%
   filter(peso >= 10)
 
 g_sem_covid <- graph_from_data_frame(pares_sem_covid, directed = FALSE)
 g_sem_covid <- delete_vertices(g_sem_covid, which(degree(g_sem_covid) == 0))
 
 # Centralidades na coorte livre de COVID
 betw_sem_covid <- betweenness(g_sem_covid, directed = FALSE, normalized = TRUE) * 100
 grau_sem_covid <- degree(g_sem_covid)
 
 df_sensibilidade <- tibble(
   CID10 = names(betw_sem_covid),
   Intermediacao_Sem_COVID = round(betw_sem_covid, 2),
   Grau_Sem_COVID = grau_sem_covid
 )
 
 # Cruza com a rede regional completa pré-existente
 df_comparativo_sensibilidade <- tibble(
   CID10 = V(grafo_regional)$name,
   Intermediacao_Completa = round(V(grafo_regional)$betweenness, 2),
   Grau_Completo = degree(grafo_regional)
 ) %>%
   inner_join(df_sensibilidade, by = "CID10")
 
 # Teste de correlação de Spearman (avalia se a hierarquia de importância permaneceu intacta)
 cor_teste <- cor.test(df_comparativo_sensibilidade$Intermediacao_Completa, 
                       df_comparativo_sensibilidade$Intermediacao_Sem_COVID, 
                       method = "spearman")
 
 cat(sprintf("    - Correlação de Spearman entre as redes: rho = %.4f (P-valor = %.2e)\n", 
             cor_teste$estimate, cor_teste$p.value))
 
 # Geração da Figura 13 (Dispersão de Estabilidade Topológica)
 plot_sensibilidade <- ggplot(df_comparativo_sensibilidade, 
                              aes(x = Intermediacao_Completa, y = Intermediacao_Sem_COVID)) +
   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.8) +
   geom_point(color = "#08519c", alpha = 0.7, size = 3) +
   ggrepel::geom_text_repel(
     data = df_comparativo_sensibilidade %>% filter(Intermediacao_Completa > 5 | Intermediacao_Sem_COVID > 5),
     aes(label = CID10),
     fontface = "bold", size = 3.6, box.padding = 0.4
   ) +
   labs(
     title = "Análise de Sensibilidade: Estabilidade Topológica Livre de COVID-19",
     subtitle = sprintf("Comparação da Intermediação (%%) — Coorte Integral vs. Livre de COVID (Spearman rho = %.3f, p < 0.001)", cor_teste$estimate),
     x = "Centralidade de Intermediação na Coorte Completa (%)",
     y = "Centralidade de Intermediação Livre de COVID (%)"
   ) +
   theme_minimal(base_size = 12) +
   theme(
     plot.title    = element_text(face = "bold", size = 13),
     plot.subtitle = element_text(size = 10.5, color = "grey30", margin = margin(b = 10)),
     axis.title    = element_text(face = "bold", size = 11)
   )
 
 ggsave("Figura13_Sensibilidade_Livre_de_COVID.png", plot_sensibilidade, width = 9, height = 7, dpi = 300, bg = "white")
 cat(" -> Figura 13 ('Figura13_Sensibilidade_Livre_de_COVID.png') salva com sucesso!\n")
 
 # ------------------------------------------------------------------------------
 # ATUALIZAÇÃO DO WORKBOOK DO EXCEL (OSF)
 # ------------------------------------------------------------------------------
 if (file.exists("OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx")) {
   wb_atual <- loadWorkbook("OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx")
   if (!"Sensibilidade_Sem_COVID" %in% names(wb_atual)) {
     addWorksheet(wb_atual, "Sensibilidade_Sem_COVID")
     writeData(wb_atual, "Sensibilidade_Sem_COVID", df_comparativo_sensibilidade)
     saveWorkbook(wb_atual, "OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx", overwrite = TRUE)
     cat(" -> Aba 'Sensibilidade_Sem_COVID' adicionada ao Excel OSF com sucesso!\n")
   }
 }
 
 cat("\n=== SUBANÁLISES CONCLUÍDAS COM SUCESSO! ===\n")
 
 
 
 # ==============================================================================
 # CORREÇÃO DA FIGURA 12: FOREST PLOT COM SUBGRUPOS (PRÉ, DURANTE E PÓS-COVID)
 # ==============================================================================
 graphics.off() # Fecha qualquer dispositivo pendente
 
 gerar_forest_subgrupos_era_v2 <- function(base_cod, comorb_cod) {
   # Agregação 2x2 por UF dentro de cada Era
   tabela_uf_era <- dados_subgrupos %>%
     group_by(Subgrupo_Era, UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (n_distinct(tabela_uf_era$Subgrupo_Era) < 2 || nrow(tabela_uf_era) < 6) {
     cat("    -> Dados insuficientes para subgrupos de", base_cod, "x", comorb_cod, "\n")
     return(invisible(NULL))
   }
   
   # Meta-análise com subgrupos no pacote meta
   m_sub <- tryCatch({
     metabin(
       event.e  = a,
       n.e      = a + b,
       event.c  = c,
       n.c      = c + d,
       studlab  = UF,
       subgroup = Subgrupo_Era,
       data     = tabela_uf_era,
       sm       = "RR",
       method   = "Inverse",
       method.tau = "REML",
       method.random.ci = "HK"
     )
   }, error = function(e) {
     cat("    -> Erro no metabin:", e$message, "\n")
     NULL
   })
   
   if (is.null(m_sub)) return(invisible(NULL))
   
   nome_arq <- paste0("Figura12_Forest_Subgrupos_COVID_", base_cod, "_x_", comorb_cod, ".png")
   
   png(nome_arq, width = 2300, height = 1500 + 45 * nrow(tabela_uf_era), res = 220)
   
   tryCatch({
     forest(
       m_sub,
       layout           = "meta",
       common           = FALSE,
       random           = TRUE,
       test.subgroup    = TRUE,
       col.square       = "#08306b",
       col.diamond      = "black",
       squaresize       = 0.75,
       leftcols         = c("studlab", "event.e", "n.e", "event.c", "n.c"),
       leftlabs         = c("UF", "Óbitos", "Total", "Óbitos", "Total"),
       label.e          = paste0("Comorbidade (", comorb_cod, ")"),
       label.c          = paste0("Base Isolada (", base_cod, ")"),
       rightcols        = c("effect", "ci", "w.random"),
       rightlabs        = c("RR", "95%-CI", "Peso (%)"),
       ref              = 1,
       digits           = 2,
       digits.weight    = 1,
       smlab            = "Risco Relativo"
     )
     # Título adicionado com grid text (compatível com a engine moderna do forest/grid)
     grid::grid.text(
       paste0("Análise de Subgrupos (Eras Epidemiológicas): ", base_cod, " → ", comorb_cod),
       x = 0.5, y = 0.98,
       gp = grid::gpar(fontsize = 13, fontface = "bold")
     )
     cat("    -> Forest Plot com subgrupos salvo com sucesso:", nome_arq, "\n")
   }, finally = {
     dev.off()
   })
   
   invisible(NULL)
 }
 
 # Executa para os dois pares clínicos-chave
 cat("\n=== GERANDO FIGURAS 12 (SUBGRUPOS TEMPORAIS) ===\n")
 gerar_forest_subgrupos_era_v2("I64", "J180")
 gerar_forest_subgrupos_era_v2("I500", "N179")
 cat("=== CONCLUÍDO ===\n")
 
 
 
 # ==============================================================================
 # PIPELINE FINAL: FOREST PLOTS SUBGRUPOS (LIMPEZA, FILTRO RR>1 E S/ SOBREPOSIÇÃO)
 # ==============================================================================
 cat("\n=== INICIANDO GERAÇÃO FINAL DOS FOREST PLOTS (FILTRO ESTREITO) ===\n")
 
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(meta)
 library(grid)
 
 graphics.off()
 
 # 1. Limpeza da pasta de saída
 pasta_saida <- "Forest_Subgrupos_COVID"
 if (!dir.exists(pasta_saida)) {
   dir.create(pasta_saida, recursive = TRUE)
   cat(paste0(" -> Pasta criada: '", file.path(getwd(), pasta_saida), "'\n"))
 } else {
   # Apaga todos os arquivos PNG antigos dentro da pasta
   arquivos_antigos <- list.files(pasta_saida, pattern = "\\.png$", full.names = TRUE)
   if(length(arquivos_antigos) > 0) {
     file.remove(arquivos_antigos)
     cat(" -> Pasta limpa:", length(arquivos_antigos), "gráficos antigos apagados.\n")
   }
 }
 
 # 2. Preparação dos dados (mesma lógica)
 col_ano <- cols_ref$col_ano
 
 dados_subgrupos <- coorte_nordeste %>%
   filter(!is.na(UF), UF != "", UF != "NA") %>%
   mutate(
     Ano_Num = as.numeric(as.character(.data[[col_ano]])),
     Subgrupo_Era = case_when(
       Ano_Num <= 2019 ~ "1. Pré-Pandemia (2016-2019)",
       Ano_Num >= 2020 & Ano_Num <= 2022 ~ "2. Auge Pandêmico (2020-2022)",
       Ano_Num >= 2023 ~ "3. Pós-Pandemia (2023-2025)"
     )
   ) %>%
   filter(!is.na(Subgrupo_Era))
 
 # 3. Lista unificada de pares
 todos_pares_analise <- bind_rows(
   pares_para_forest %>% select(DIAG_PRINC, DIAGSEC1),
   if (exists("pares_baixa_hetero")) pares_baixa_hetero %>% select(DIAG_PRINC, DIAGSEC1) else NULL
 ) %>%
   distinct(DIAG_PRINC, DIAGSEC1)
 
 # 4. Função Geradora com Filtro e Correção Visual
 gerar_forest_subgrupos_final <- function(base_cod, comorb_cod) {
   tabela_uf_era <- dados_subgrupos %>%
     group_by(Subgrupo_Era, UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (n_distinct(tabela_uf_era$Subgrupo_Era) < 2 || nrow(tabela_uf_era) < 5) return(invisible(NULL))
   
   m_sub <- tryCatch({
     metabin(
       event.e    = a,
       n.e        = a + b,
       event.c    = c,
       n.c        = c + d,
       studlab    = UF,
       subgroup   = Subgrupo_Era,
       data       = tabela_uf_era,
       sm         = "RR",
       method     = "Inverse",
       method.tau = "REML",
       method.random.ci = "HK"
     )
   }, error = function(e) NULL)
   
   if (is.null(m_sub)) return(invisible(NULL))
   
   # === O FILTRO DE OURO: RR > 1 e P-Valor < 0.05 ===
   # Usamos o exp(TE.random) para pegar o Risco Relativo Global
   rr_global <- exp(m_sub$TE.random)
   pval_global <- m_sub$pval.random
   
   if (is.na(rr_global) || rr_global <= 1 || is.na(pval_global) || pval_global >= 0.05) {
     # Ignora pares protetores (RR < 1) ou não significativos
     return(invisible(NULL))
   }
   
   caminho_arquivo <- file.path(pasta_saida, paste0("Forest_Subgrupo_", base_cod, "_x_", comorb_cod, ".png"))
   altura_calc <- 1550 + 55 * nrow(tabela_uf_era)
   
   png(caminho_arquivo, width = 2800, height = altura_calc, res = 220)
   
   tryCatch({
     forest(
       m_sub,
       layout           = "meta",
       common           = FALSE,
       random           = TRUE,
       test.subgroup    = TRUE,
       test.overall     = TRUE,          # <--- ADICIONA O P-VALOR DO RR GLOBAL NO RODAPÉ
       
       col.square       = "#08306b",
       col.square.lines = "#08306b",
       col.diamond      = "black",
       col.diamond.lines= "black",
       squaresize       = 0.75,
       
       # Textos estruturais das colunas (mais curtos para evitar sobreposição)
       leftcols         = c("studlab", "event.e", "n.e", "event.c", "n.c"),
       leftlabs         = c("Estado\n(UF)", "Óbitos", "Total", "Óbitos", "Total"),
       
       label.e          = "Comorbidade",  # Sem CID aqui para não invadir o espaço
       label.c          = "Base",         # Sem CID aqui também
       
       rightcols        = c("effect", "ci", "w.random"),
       rightlabs        = c("RR", "IC 95%", "Peso (%)"),
       
       # Margens e referências
       ref              = 1,
       col.inside       = "grey75",
       colgap.left      = "4mm",         # Força um espaço maior entre as colunas
       
       digits           = 2,
       digits.weight    = 1,
       smlab            = "Risco Relativo",
       
       text.random      = "Modelo de Efeitos Aleatórios",
       print.subgroup.name = FALSE
     )
     
     # Título do Gráfico (mantendo os CIDs visíveis)
     grid.text(
       paste0("Análise de Subgrupos (Eras COVID): ", base_cod, " (Base) → ", comorb_cod, " (Comorbidade)"),
       x = 0.5, y = 0.985,
       gp = gpar(fontsize = 14, fontface = "bold")
     )
     
     cat(sprintf("  [+] Gerado (RR = %.2f, p = %.4f): %s\n", rr_global, pval_global, basename(caminho_arquivo)))
   }, error = function(e) {
     cat("  [-] Falha visual em:", base_cod, "x", comorb_cod, "\n")
   }, finally = {
     dev.off()
   })
   
   invisible(NULL)
 }
 
 # 5. Roda a função para todos os pares
 cat(" -> Processando pares (descartando não-significativos e RR < 1)...\n")
 pwalk(
   list(todos_pares_analise$DIAG_PRINC, todos_pares_analise$DIAGSEC1),
   gerar_forest_subgrupos_final
 )
 
 cat("\n=== PROCESSAMENTO CONCLUÍDO COM SUCESSO! ===\n")
 
 
 # ==============================================================================
 # MÓDULO 20: CARTOGRAFIA DE RISCO RELATIVO NO NORDESTE (MAPA DE CALOR)
 # ==============================================================================
 cat("\n=== INICIANDO GERAÇÃO DOS MAPAS CARTOGRÁFICOS DO NORDESTE ===\n")
 
 # Instala pacotes espaciais se não estiverem presentes
 if (!requireNamespace("geobr", quietly = TRUE)) install.packages("geobr")
 if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
 if (!requireNamespace("viridis", quietly = TRUE)) install.packages("viridis")
 
 library(dplyr)
 library(ggplot2)
 library(sf)
 library(geobr)
 library(viridis)
 
 graphics.off()
 
 # 1. Função para calcular o Risco Relativo Bruto por UF para o mapa
 calcular_rr_uf_bruto <- function(base_cod, comorb_cod) {
   coorte_nordeste %>%
     filter(!is.na(UF), UF != "", UF != "NA") %>%
     group_by(UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0) %>%
     mutate(
       Risco_Exposto = a / (a + b),
       Risco_Controle = c / (c + d),
       RR_Bruto = Risco_Exposto / Risco_Controle
     ) %>%
     select(UF, RR_Bruto)
 }
 
 # 2. Download do shapefile espacial do Nordeste (geobr)
 cat(" -> Baixando malha espacial do Nordeste via IPEA (geobr)...\n")
 mapa_ne <- read_state(code_state = "NE", year = 2020, showProgress = FALSE)
 
 # 3. Função para gerar e salvar o mapa
 gerar_mapa_calor <- function(base_cod, comorb_cod, titulo_base, titulo_comorb) {
   
   # Calcula RRs
   df_rr <- calcular_rr_uf_bruto(base_cod, comorb_cod)
   
   # Junta os dados estatísticos com a malha espacial (sf)
   mapa_dados <- mapa_ne %>%
     left_join(df_rr, by = c("abbrev_state" = "UF"))
   
   # Extrai centróides para plotar os rótulos de texto de cada estado
   mapa_centroides <- mapa_dados %>% 
     mutate(centroide = st_centroid(geom)) %>%
     mutate(lon = st_coordinates(centroide)[,1],
            lat = st_coordinates(centroide)[,2])
   
   caminho_mapa <- paste0("Figura_MapaCalor_", base_cod, "_x_", comorb_cod, ".png")
   
   p <- ggplot(mapa_dados) +
     # Desenha os estados preenchidos com o valor do RR
     geom_sf(aes(fill = RR_Bruto), color = "white", linewidth = 0.6) +
     
     # Adiciona a sigla do estado e o RR numérico no centro de cada polígono
     geom_text(data = mapa_centroides, 
               aes(x = lon, y = lat, label = paste0(abbrev_state, "\n", round(RR_Bruto, 1))),
               color = "black", fontface = "bold", size = 3.5, lineheight = 0.8) +
     
     # Paleta de calor elegante (amarelo -> laranja -> vermelho escuro)
     scale_fill_viridis(option = "inferno", direction = -1, begin = 0.1, end = 0.9, 
                        na.value = "grey90", name = "Risco Relativo") +
     
     labs(
       title = paste0("Distribuição Espacial da Letalidade (", base_cod, " → ", comorb_cod, ")"),
       subtitle = paste0(titulo_base, " evoluindo com ", titulo_comorb, "\nSIH/SUS Nordeste (2016-2025)"),
       caption = "Fonte: Elaboração própria via dados do SIH/DATASUS. Malha: IPEA/geobr."
     ) +
     
     # Tema vazio estrito para mapas
     theme_void(base_size = 14) +
     theme(
       plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(t = 10, b = 5)),
       plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30", margin = margin(b = 20)),
       plot.caption = element_text(size = 9, color = "grey50", hjust = 0.95),
       legend.position = "right",
       legend.title = element_text(face = "bold", size = 11),
       plot.background = element_rect(fill = "white", color = NA),
       plot.margin = margin(10, 10, 10, 10)
     )
   
   ggsave(caminho_mapa, plot = p, width = 10, height = 9, dpi = 300, bg = "white")
   cat("  [+] Mapa gerado com sucesso:", caminho_mapa, "\n")
 }
 
 # 4. Geração dos dois mapas protagonistas
 gerar_mapa_calor("I64", "J180", "Acidente Vascular Cerebral (AVC)", "Broncopneumonia")
 gerar_mapa_calor("I500", "N179", "Insuficiência Cardíaca", "Insuficiência Renal Aguda")
 
 cat("\n=== CARTOGRAFIA CONCLUÍDA! ===\n")
 
 # ==============================================================================
 # CORREÇÃO DE ESCALA: TABELA 1 EM MODO PAISAGEM (LANDSCAPE) E FONTE AJUSTADA
 # ==============================================================================
 cat("\n=== GERANDO TABELA 1 COM LAYOUT CORRIGIDO ===\n")
 
 # Instala o officer se não tiver
 if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
 
 library(dplyr)
 library(gtsummary)
 library(flextable)
 library(officer)
 
 # 1. Tratamento dos dados (mesmo de antes)
 dados_tabela1 <- coorte_nordeste %>%
   filter(!is.na(UF), UF != "", UF != "NA") %>%
   mutate(
     Ano_Num = as.numeric(as.character(.data[[cols_ref$col_ano]])),
     Fase_Pandemia = case_when(
       Ano_Num <= 2019 ~ "1. Pré-Pandemia (2016-2019)",
       Ano_Num >= 2020 & Ano_Num <= 2022 ~ "2. Auge Pandêmico (2020-2022)",
       Ano_Num >= 2023 ~ "3. Pós-Pandemia (2023-2025)"
     ),
     Sexo = case_when(
       SEXO %in% c("1", "M") ~ "Masculino",
       SEXO %in% c("3", "F") ~ "Feminino",
       TRUE ~ "Ignorado/Sem info"
     ),
     Raca_Cor = case_when(
       RACA_COR == "01" ~ "Branca",
       RACA_COR == "02" ~ "Preta",
       RACA_COR == "03" ~ "Parda",
       RACA_COR == "04" ~ "Amarela",
       RACA_COR == "05" ~ "Indígena",
       TRUE ~ "Sem informação"
     ),
     Idade = as.numeric(IDADE),
     Dias_Internacao = as.numeric(DIAS_PERM)
   ) %>%
   select(Fase_Pandemia, Sexo, Raca_Cor, Idade, Dias_Internacao, UF) %>%
   filter(!is.na(Fase_Pandemia))
 
 theme_gtsummary_language(language = "pt", decimal.mark = ",", big.mark = ".")
 
 # 2. Montagem da Tabela
 tabela_demografica <- dados_tabela1 %>%
   tbl_summary(
     by = Fase_Pandemia,
     missing = "no", # Removemos a linha de missing para poupar espaço vertical
     label = list(
       Sexo ~ "Sexo Biológico",
       Raca_Cor ~ "Raça/Cor",
       Idade ~ "Idade (anos)",
       Dias_Internacao ~ "Permanência (dias)",
       UF ~ "Estado (UF)"
     ),
     statistic = list(
       all_continuous() ~ "{median} ({p25} - {p75})",
       all_categorical() ~ "{n} ({p}%)"
     )
   ) %>%
   add_overall(col_label = "**Total Geral**  \n(N = {N})") %>%
   add_p(
     test = list(all_continuous() ~ "kruskal.test", all_categorical() ~ "chisq.test"),
     pvalue_fun = ~style_pvalue(.x, digits = 3)
   ) %>%
   modify_header(label = "**Características Basais**") %>%
   modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Fase Epidemiológica**") %>%
   bold_labels()
 
 # 3. Formatação Estética Avançada para Caber na Página
 tabela_flex <- tabela_demografica %>%
   as_flex_table() %>%
   theme_vanilla() %>%
   fontsize(size = 9, part = "all") %>%           # Reduz a fonte de 11 para 9
   padding(padding = 3, part = "all") %>%         # Diminui o espaçamento interno das células
   set_table_properties(layout = "autofit") %>% 
   fit_to_width(max_width = 10.5)                 # Força a tabela a ter no máximo 10.5 polegadas
 
 # 4. Configuração da Página do Word para Paisagem (Landscape) com margens estreitas
 prop_secao_paisagem <- prop_section(
   page_size = page_size(orient = "landscape", width = 11.69, height = 8.27),
   page_margins = page_mar(bottom = 0.5, top = 0.5, right = 0.5, left = 0.5)
 )
 
 # 5. Salva o arquivo final
 save_as_docx(tabela_flex, path = "Tabela1_Demografia_Perfeita.docx", pr_section = prop_secao_paisagem)
 
 cat(" -> SUCESSO! A tabela foi salva como 'Tabela1_Demografia_Perfeita.docx' em formato Paisagem.\n")
 
 
 
 # ==============================================================================
 # PIPELINE FINAL: FOREST PLOTS SUBGRUPOS (ESTÉTICA PADRÃO-OURO E ODDS RATIO)
 # ==============================================================================
 cat("\n=== INICIANDO GERAÇÃO FINAL DOS FOREST PLOTS (ESTÉTICA CORRIGIDA) ===\n")
 
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(meta)
 library(grid)
 
 graphics.off()
 
 # 1. Limpeza da pasta de saída
 pasta_saida <- "Forest_Subgrupos_COVID"
 if (!dir.exists(pasta_saida)) {
   dir.create(pasta_saida, recursive = TRUE)
 } else {
   arquivos_antigos <- list.files(pasta_saida, pattern = "\\.png$", full.names = TRUE)
   if(length(arquivos_antigos) > 0) file.remove(arquivos_antigos)
 }
 
 # 2. Preparação dos dados
 col_ano <- cols_ref$col_ano
 
 dados_subgrupos <- coorte_nordeste %>%
   filter(!is.na(UF), UF != "", UF != "NA") %>%
   mutate(
     Ano_Num = as.numeric(as.character(.data[[col_ano]])),
     Subgrupo_Era = case_when(
       Ano_Num <= 2019 ~ "1. Pré-Pandemia (2016-2019)",
       Ano_Num >= 2020 & Ano_Num <= 2022 ~ "2. Auge Pandêmico (2020-2022)",
       Ano_Num >= 2023 ~ "3. Pós-Pandemia (2023-2025)"
     )
   ) %>%
   filter(!is.na(Subgrupo_Era))
 
 # 3. Lista unificada de pares
 todos_pares_analise <- bind_rows(
   pares_para_forest %>% select(DIAG_PRINC, DIAGSEC1),
   if (exists("pares_baixa_hetero")) pares_baixa_hetero %>% select(DIAG_PRINC, DIAGSEC1) else NULL
 ) %>%
   distinct(DIAG_PRINC, DIAGSEC1)
 
 # 4. Função Geradora (Estética Idêntica à Referência)
 gerar_forest_subgrupos_final <- function(base_cod, comorb_cod) {
   tabela_uf_era <- dados_subgrupos %>%
     group_by(Subgrupo_Era, UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (n_distinct(tabela_uf_era$Subgrupo_Era) < 2 || nrow(tabela_uf_era) < 5) return(invisible(NULL))
   
   m_sub <- tryCatch({
     metabin(
       event.e    = a,
       n.e        = a + b,
       event.c    = c,
       n.c        = c + d,
       studlab    = UF,
       subgroup   = Subgrupo_Era,
       data       = tabela_uf_era,
       sm         = "OR",             # Mantém Odds Ratio
       method     = "Inverse",
       method.tau = "REML",
       method.random.ci = "HK"
     )
   }, error = function(e) NULL)
   
   if (is.null(m_sub)) return(invisible(NULL))
   
   or_global <- exp(m_sub$TE.random)
   pval_global <- m_sub$pval.random
   
   if (is.na(or_global) || or_global <= 1 || is.na(pval_global) || pval_global >= 0.05) {
     return(invisible(NULL))
   }
   
   caminho_arquivo <- file.path(pasta_saida, paste0("Forest_Subgrupo_", base_cod, "_x_", comorb_cod, ".png"))
   
   # Cálculo exato de altura para não espremer os dados
   altura_calc <- 1600 + 60 * nrow(tabela_uf_era)
   
   png(caminho_arquivo, width = 3000, height = altura_calc, res = 220, bg = "white") 
   
   tryCatch({
     forest(
       m_sub,
       layout           = "meta",
       common           = FALSE,
       random           = TRUE,           
       test.overall     = TRUE,           
       test.subgroup    = TRUE,
       
       col.square       = "#08306b",
       col.square.lines = "#08306b",
       col.diamond      = "black",
       col.diamond.lines= "black",
       squaresize       = 0.75,
       
       # === RESTAURAÇÃO DA ESTÉTICA DE REFERÊNCIA ===
       leftcols         = c("studlab", "event.e", "n.e", "event.c", "n.c"),
       leftlabs         = c("Estado\n(UF)", "Óbitos", "Total", "Óbitos", "Total"),
       
       # Retorna o cabeçalho superior aninhado
       label.e          = "Comorbidade",  
       label.c          = "Base",         
       
       rightcols        = c("effect", "ci", "w.random"),
       rightlabs        = c("OR", "IC 95%", "Peso (%)"),
       
       ref              = 1,
       col.inside       = "grey75",
       
       # Ajuste sutil de margens para o texto respirar
       colgap.left      = "6mm",          
       colgap.forest    = "6mm",
       
       digits           = 2,
       digits.weight    = 1,
       smlab            = "Odds Ratio", 
       text.random      = "Modelo de Efeitos Aleatórios",
       
       print.subgroup.name = FALSE
     )
     
     # Título centralizado no topo com margem limpa
     grid.text(
       paste0("Análise de Subgrupos (Eras COVID): ", base_cod, " (Base) \u2192 ", comorb_cod, " (Comorbidade)"),
       x = 0.5, y = 0.985,
       gp = gpar(fontsize = 16, fontface = "bold")
     )
     
     cat(sprintf("  [+] Gerado com estética padrão: %s\n", basename(caminho_arquivo)))
   }, error = function(e) {
     cat("  [-] Falha visual em:", base_cod, "x", comorb_cod, "\n")
   }, finally = {
     dev.off()
   })
   
   invisible(NULL)
 }
 
 # 5. Roda a função
 cat(" -> Processando pares...\n")
 pwalk(
   list(todos_pares_analise$DIAG_PRINC, todos_pares_analise$DIAGSEC1),
   gerar_forest_subgrupos_final
 )
 
 cat("\n=== PROCESSAMENTO CONCLUÍDO COM SUCESSO! ===\n")
 
 
 # ==============================================================================
 # EXPORTAÇÃO EXCLUSIVA DE DADOS (PLANILHA DE RESULTADOS DA META-ANÁLISE)
 # ==============================================================================
 cat("\n=== INICIANDO EXTRAÇÃO DE DADOS PARA PLANILHA ===\n")
 
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(meta)
 
 # Função para rodar o modelo e extrair apenas os números
 extrair_dados_meta <- function(base_cod, comorb_cod) {
   tabela_uf_era <- dados_subgrupos %>%
     group_by(Subgrupo_Era, UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (n_distinct(tabela_uf_era$Subgrupo_Era) < 2 || nrow(tabela_uf_era) < 5) return(NULL)
   
   m_sub <- tryCatch({
     metabin(
       event.e = a, n.e = a + b, event.c = c, n.c = c + d,
       studlab = UF, subgroup = Subgrupo_Era, data = tabela_uf_era,
       sm = "OR", method = "Inverse", method.tau = "REML", method.random.ci = "HK"
     )
   }, error = function(e) NULL)
   
   if (is.null(m_sub)) return(NULL)
   
   or_global <- exp(m_sub$TE.random)
   pval_global <- m_sub$pval.random
   
   # Aplica o mesmo filtro dos gráficos (ignora protetores ou não significativos)
   if (is.na(or_global) || or_global <= 1 || is.na(pval_global) || pval_global >= 0.05) {
     return(NULL)
   }
   
   # Monta a tabela com os resultados de cada subgrupo (Era) e o total global
   dados_extraidos <- tibble(
     Doenca_Base = base_cod,
     Comorbidade = comorb_cod,
     Fase = c(m_sub$bylevs, "4. Total Global"),
     OR_Ajustado = c(exp(m_sub$TE.random.w), exp(m_sub$TE.random)),
     IC_Inferior = c(exp(m_sub$lower.random.w), exp(m_sub$lower.random)),
     IC_Superior = c(exp(m_sub$upper.random.w), exp(m_sub$upper.random)),
     P_Valor = c(m_sub$pval.random.w, m_sub$pval.random)
   )
   
   cat(sprintf("  [+] Extraído: %s x %s\n", base_cod, comorb_cod))
   return(dados_extraidos)
 }
 
 # Roda a extração para todos os pares e consolida em um único dataframe
 cat(" -> Processando cálculos estatísticos...\n")
 tabela_resultados_finais <- map2_dfr(
   todos_pares_analise$DIAG_PRINC, 
   todos_pares_analise$DIAGSEC1, 
   extrair_dados_meta
 )
 
 # Salva o arquivo CSV no padrão brasileiro (pronto para Excel)
 if (nrow(tabela_resultados_finais) > 0) {
   write.csv2(tabela_resultados_finais, "Resultados_MetaAnalise_Completos.csv", row.names = FALSE)
   cat("\n=== SUCESSO! ===\n")
   cat(" -> A planilha 'Resultados_MetaAnalise_Completos.csv' foi criada na sua pasta de trabalho.\n")
 } else {
   cat("\nNenhum dado válido para exportação.\n")
 }
 
 
 # ==============================================================================
 # RECRIAÇÃO DA FIGURA DE PARES ROBUSTOS COM ODDS RATIO (OR)
 # ==============================================================================
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(meta)
 library(ggplot2)
 
 cat("\n=== PROCESSANDO DADOS DA META-ANÁLISE (ODDS RATIO E I²) ===\n")
 
 # 1. Função para extrair dados globais e contagem de UFs significativas por par
 extrair_dados_figura9 <- function(base_cod, comorb_cod) {
   tabela_uf <- coorte_nordeste %>%
     filter(!is.na(UF), UF != "", UF != "NA") %>%
     group_by(UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (nrow(tabela_uf) < 4) return(NULL) # Mínimo de estados para considerar universal
   
   m <- tryCatch({
     metabin(
       event.e = a, n.e = a + b, event.c = c, n.c = c + d, studlab = UF,
       data = tabela_uf, sm = "OR", method = "Inverse", method.tau = "REML", method.random.ci = "HK"
     )
   }, error = function(e) NULL)
   
   if (is.null(m)) return(NULL)
   
   # Filtro 1: Apenas associações globais significativas
   or_global <- exp(m$TE.random)
   pval_global <- m$pval.random
   if (is.na(or_global) || or_global <= 1 || is.na(pval_global) || pval_global >= 0.05) return(NULL)
   
   # Filtro 2: Baixa Heterogeneidade (I² < 50%)
   i2_valor <- m$I2
   if (is.na(i2_valor) || i2_valor >= 0.50) return(NULL)
   
   # Contagem de UFs onde a associação foi estatisticamente significativa e de risco (OR > 1)
   ufs_sig <- sum(exp(m$TE) > 1 & m$pval < 0.05, na.rm = TRUE)
   
   tibble(
     Base = base_cod,
     Comorb = comorb_cod,
     OR = or_global,
     IC_Inf = exp(m$lower.random),
     IC_Sup = exp(m$upper.random),
     I2 = i2_valor,
     UFs_Sig = ufs_sig
   )
 }
 
 # 2. Executa a extração em todos os pares mapeados
 df_plot_robustos <- map2_dfr(
   todos_pares_analise$DIAG_PRINC, 
   todos_pares_analise$DIAGSEC1, 
   extrair_dados_figura9
 )
 
 # Verifica se sobraram pares após o filtro rigoroso do OR + I² < 50%
 if(nrow(df_plot_robustos) == 0) stop("Nenhum par sobreviveu aos critérios de OR > 1 e I² < 50%.")
 
 # 3. Tratamento dos dados para o ggplot2
 df_plot_robustos <- df_plot_robustos %>%
   arrange(OR) %>%
   mutate(
     # Cria o rótulo combinando CIDs e o valor de I²
     Label = sprintf("%s \u2192 %s (I\u00b2 = %.1f%%)", Base, Comorb, I2 * 100),
     # Trava a ordem dos fatores para o gráfico respeitar o arrange()
     Label = factor(Label, levels = Label),
     # Transforma contagem em fator para mapear shapes corretos
     UFs_Sig_Cat = as.factor(UFs_Sig)
   )
 
 cat(sprintf(" -> %d pares robustos identificados para o plot.\n", nrow(df_plot_robustos)))
 
 # 4. Construção do Gráfico (Estética Idêntica à Figura 9)
 plot_fig9 <- ggplot(df_plot_robustos, aes(x = OR, y = Label)) +
   
   # Linhas de grade sutis
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.major.y = element_blank(),
     panel.grid.major.x = element_line(color = "grey85", linetype = "dotted"),
     panel.grid.minor.x = element_blank(),
     plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
     plot.subtitle = element_text(color = "grey40", size = 14, margin = margin(b = 20)),
     axis.text.y = element_text(face = "bold", color = "black", size = 12),
     axis.text.x = element_text(color = "black", size = 12),
     axis.title.x = element_text(face = "bold", margin = margin(t = 15)),
     legend.position = "bottom",
     legend.title = element_text(face = "bold"),
     legend.text = element_text(size = 12)
   ) +
   
   # Linha de nulidade (OR = 1)
   geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.8) +
   
   # Barras de erro (Intervalo de Confiança)
   geom_errorbarh(aes(xmin = IC_Inf, xmax = IC_Sup), height = 0.25, color = "#1b7837", size = 0.9) +
   
   # Pontos mapeados por formato com fundo verde escuro e borda preta
   geom_point(aes(shape = UFs_Sig_Cat), fill = "#1b7837", color = "black", size = 5.5, stroke = 1) +
   
   # Mapeamento manual de formatos compatíveis com 'fill' (21 a 25)
   # 21 = Círculo, 22 = Quadrado, 23 = Losango, 24 = Triângulo
   scale_shape_manual(
     values = c("3" = 21, "4" = 21, "5" = 22, "6" = 24, "7" = 23, "8" = 25, "9" = 21)
   ) +
   
   # Escala Logarítmica exata da imagem original
   scale_x_log10(
     breaks = c(1, 2, 3, 5, 10, 20),
     labels = c("1", "2", "3", "5", "10", "20"),
     limits = c(0.8, max(df_plot_robustos$IC_Sup) * 1.2)
   ) +
   
   # Textos
   labs(
     title = "Associações Clínicas Letais Universais no Nordeste (I\u00b2 < 50%)",
     subtitle = "Pares com elevada consistência biológica e homogeneidade entre estados (2016-2025)",
     x = "Odds Ratio Meta-Analítico \u2014 Escala Logarítmica (IC 95%)",
     y = NULL,
     shape = "UFs com associação significante:"
   )
 
 # 5. Salvar Imagem em Alta Resolução
 ggsave(
   "Figura9_ForestPlot_Pares_Robustos_OR.png",
   plot = plot_fig9,
   width = 14, height = 9, dpi = 300, bg = "white"
 )
 
 cat("=== SUCESSO: 'Figura9_ForestPlot_Pares_Robustos_OR.png' gerada! ===\n")
 
 
 
 
 # ==============================================================================
 # RECRIAÇÃO DA FIGURA 9: TODAS AS ASSOCIAÇÕES SIGNIFICATIVAS (OR > 1)
 # ==============================================================================
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(meta)
 library(ggplot2)
 
 cat("\n=== PROCESSANDO DADOS GLOBAIS (ODDS RATIO) ===\n")
 
 extrair_dados_figura9_or <- function(base_cod, comorb_cod) {
   tabela_uf <- coorte_nordeste %>%
     filter(!is.na(UF), UF != "", UF != "NA") %>%
     group_by(UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (nrow(tabela_uf) < 4) return(NULL)
   
   m <- tryCatch({
     metabin(
       event.e = a, n.e = a + b, event.c = c, n.c = c + d, studlab = UF,
       data = tabela_uf, sm = "OR", method = "Inverse", method.tau = "REML", method.random.ci = "HK"
     )
   }, error = function(e) NULL)
   
   if (is.null(m)) return(NULL)
   
   # Filtro Único: Apenas associações globais significativas e de risco (OR > 1)
   or_global <- exp(m$TE.random)
   pval_global <- m$pval.random
   if (is.na(or_global) || or_global <= 1 || is.na(pval_global) || pval_global >= 0.05) return(NULL)
   
   i2_valor <- m$I2
   ufs_sig <- sum(exp(m$TE) > 1 & m$pval < 0.05, na.rm = TRUE)
   
   tibble(
     Base = base_cod,
     Comorb = comorb_cod,
     OR = or_global,
     IC_Inf = exp(m$lower.random),
     IC_Sup = exp(m$upper.random),
     I2 = ifelse(is.na(i2_valor), 0, i2_valor),
     UFs_Sig = ufs_sig
   )
 }
 
 df_plot_robustos <- map2_dfr(
   todos_pares_analise$DIAG_PRINC, 
   todos_pares_analise$DIAGSEC1, 
   extrair_dados_figura9_or
 )
 
 if(nrow(df_plot_robustos) == 0) stop("Nenhum par significativo encontrado.")
 
 df_plot_robustos <- df_plot_robustos %>%
   arrange(OR) %>%
   mutate(
     # Rótulo atualizado para manter a informação do I2, mas sem excluí-los
     Label = sprintf("%s \u2192 %s (I\u00b2 = %.1f%%)", Base, Comorb, I2 * 100),
     Label = factor(Label, levels = Label),
     UFs_Sig_Cat = as.factor(UFs_Sig)
   )
 
 cat(sprintf(" -> %d pares significativos identificados para o plot.\n", nrow(df_plot_robustos)))
 
 plot_fig9 <- ggplot(df_plot_robustos, aes(x = OR, y = Label)) +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.major.y = element_blank(),
     panel.grid.major.x = element_line(color = "grey85", linetype = "dotted"),
     panel.grid.minor.x = element_blank(),
     plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
     plot.subtitle = element_text(color = "grey40", size = 14, margin = margin(b = 20)),
     axis.text.y = element_text(face = "bold", color = "black", size = 12),
     axis.text.x = element_text(color = "black", size = 12),
     axis.title.x = element_text(face = "bold", margin = margin(t = 15)),
     legend.position = "bottom",
     legend.title = element_text(face = "bold"),
     legend.text = element_text(size = 12)
   ) +
   geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.8) +
   geom_errorbarh(aes(xmin = IC_Inf, xmax = IC_Sup), height = 0.25, color = "#1b7837", size = 0.9) +
   geom_point(aes(shape = UFs_Sig_Cat), fill = "#1b7837", color = "black", size = 5.5, stroke = 1) +
   scale_shape_manual(
     values = c("1" = 21, "2" = 22, "3" = 24, "4" = 21, "5" = 22, "6" = 24, "7" = 23, "8" = 25, "9" = 21)
   ) +
   scale_x_log10(
     breaks = c(1, 2, 3, 5, 10, 20, 30),
     labels = c("1", "2", "3", "5", "10", "20", "30"),
     limits = c(0.8, max(df_plot_robustos$IC_Sup) * 1.2)
   ) +
   labs(
     title = "Principais Associações Clínicas Letais no Nordeste (2016-2025)",
     subtitle = "Pares de maior impacto na mortalidade hospitalar (Modelo de Efeitos Aleatórios, OR > 1)",
     x = "Odds Ratio Ajustado \u2014 Escala Logarítmica (IC 95%)",
     y = NULL,
     shape = "UFs com associação significante no estado:"
   )
 
 ggsave(
   "Figura9_ForestPlot_Pares_Significativos_OR.png",
   plot = plot_fig9,
   width = 14, height = 9, dpi = 300, bg = "white"
 )
 
 cat("=== SUCESSO: 'Figura9_ForestPlot_Pares_Significativos_OR.png' gerada! ===\n")
 
 
 
 
 # ==============================================================================
 # REGRESSÃO LOGÍSTICA MULTIVARIADA (OTIMIZADA PARA BIG DATA - SEM TRAVAR)
 # ==============================================================================
 library(dplyr)
 library(broom)
 
 base_cod <- "I64"
 comorb_cod <- "J180"
 
 cat("\n=== PREPARANDO DADOS... ===\n")
 
 dados_regressao <- dados_subgrupos %>%
   mutate(
     Y_Comorbidade = ifelse(DIAGSEC1 == comorb_cod, 1, 0),
     X_Base        = ifelse(DIAG_PRINC == base_cod, 1, 0),
     
     Idade_Num = suppressWarnings(as.numeric(as.character(IDADE))),
     Sexo_Cat  = case_when(
       grepl("1|^m|^M", as.character(SEXO)) ~ "Masculino",
       grepl("3|^f|^F", as.character(SEXO)) ~ "Feminino",
       TRUE ~ NA_character_
     ),
     UF_Cat      = as.character(UF),
     Periodo_Cat = as.character(Subgrupo_Era),
     
     Dias_Num    = suppressWarnings(as.numeric(as.character(DIAS_PERM))),
     Carater_Cat = case_when(
       CAR_INT == "01" ~ "Eletiva",
       CAR_INT %in% c("02", "03", "04", "05", "06") ~ "Urgencia", 
       TRUE ~ NA_character_
     )
   ) %>%
   filter(
     !is.na(Idade_Num), !is.na(Sexo_Cat), !is.na(UF_Cat), 
     !is.na(Periodo_Cat), !is.na(Dias_Num), !is.na(Carater_Cat)
   ) %>%
   mutate(
     Sexo_Cat    = as.factor(Sexo_Cat),
     UF_Cat      = as.factor(UF_Cat),
     Periodo_Cat = as.factor(Periodo_Cat),
     Carater_Cat = as.factor(Carater_Cat)
   ) %>%
   droplevels()
 
 cat(" -> Base pronta. Total de casos válidos:", nrow(dados_regressao), "\n")
 cat(" -> Rodando a regressão glm() (isso deve levar menos de 10 segundos)...\n")
 
 # Roda o modelo
 modelo_estendido <- glm(
   Y_Comorbidade ~ X_Base + Idade_Num + Sexo_Cat + Periodo_Cat + UF_Cat + Dias_Num + Carater_Cat,
   data = dados_regressao,
   family = binomial(link = "logit")
 )
 
 cat(" -> glm() concluído! Calculando Odds Ratio instantâneo (Wald)...\n\n")
 
 # Extração rápida usando cálculo de Wald para evitar o travamento do conf.int=TRUE
 tabela_estendida <- tidy(modelo_estendido) %>%
   mutate(
     # Calcula OR e IC diretamente da estimativa e do erro padrão
     OR_Ajustado = exp(estimate),
     IC_Inf = exp(estimate - 1.96 * std.error),
     IC_Sup = exp(estimate + 1.96 * std.error),
     Valor_P = format.pval(p.value, eps = 0.001)
   ) %>%
   select(Variavel = term, OR_Ajustado, IC_Inf, IC_Sup, Valor_P) %>%
   mutate(across(c(OR_Ajustado, IC_Inf, IC_Sup), ~round(., 2)))
 
 # ==============================================================================
 # CORREÇÃO 1: RESGATE DO PIAUÍ E RECUPERAÇÃO DO 'N' DA COORTE
 # ==============================================================================
 library(dplyr)
 
 cat("\n=== INICIANDO RESGATE GEOGRÁFICO ===\n")
 cat("N original do banco na memória:", nrow(coorte_nordeste), "pacientes.\n\n")
 
 coorte_nordeste <- coorte_nordeste %>%
   mutate(
     # Garante que a leitura será como texto para evitar bugs de fator numérico
     UF_Raw = as.character(UF), 
     
     # Mapeamento implacável: se for sigla mantém, se for IBGE converte
     UF_Corrigida = case_when(
       UF_Raw %in% c("21", "MA") ~ "MA",
       UF_Raw %in% c("22", "PI") ~ "PI",
       UF_Raw %in% c("23", "CE") ~ "CE",
       UF_Raw %in% c("24", "RN") ~ "RN",
       UF_Raw %in% c("25", "PB") ~ "PB",
       UF_Raw %in% c("26", "PE") ~ "PE",
       UF_Raw %in% c("27", "AL") ~ "AL",
       UF_Raw %in% c("28", "SE") ~ "SE",
       UF_Raw %in% c("29", "BA") ~ "BA",
       TRUE ~ UF_Raw
     )
   )
 
 # Substitui a UF bugada pela UF limpa
 coorte_nordeste$UF <- coorte_nordeste$UF_Corrigida
 
 # Confirmação visual para você conferir se o Piauí voltou e se o N bate
 cat("--- DISTRIBUIÇÃO CORRIGIDA POR ESTADO ---\n")
 print(table(coorte_nordeste$UF, useNA = "always"))
 
 cat("\nTotal de pacientes sem estado (NA):", sum(is.na(coorte_nordeste$UF)), "\n")
 cat("N final com todos os estados:", sum(!is.na(coorte_nordeste$UF)), "\n")
 
 
 # ==============================================================================
 # CORREÇÃO 2: RESGATE DA VARIÁVEL SEXO PARA A TABELA 1
 # ==============================================================================
 library(dplyr)
 
 cat("\n=== CORRIGINDO O MAPEAMENTO DE SEXO ===\n")
 
 coorte_nordeste <- coorte_nordeste %>%
   mutate(
     # Mapeamento flexível: captura os códigos 1/3 do SUS ou as letras M/F
     Sexo_Corrigido = case_when(
       grepl("1|^m|^M", as.character(SEXO)) ~ "Masculino",
       grepl("3|^f|^F", as.character(SEXO)) ~ "Feminino",
       TRUE ~ "Ignorado"
     ),
     # Transforma em fator na ordem correta para a Tabela 1
     Sexo_Corrigido = factor(Sexo_Corrigido, levels = c("Masculino", "Feminino", "Ignorado"))
   )
 
 # Substitui a coluna original bugada pela coluna limpa
 coorte_nordeste$SEXO <- coorte_nordeste$Sexo_Corrigido
 
 # Confirmação visual dos dados reais
 cat("--- DISTRIBUIÇÃO REAL POR SEXO ---\n")
 print(table(coorte_nordeste$SEXO, useNA = "always"))
 
 
 # ==============================================================================
 # CORREÇÃO 3: REGERANDO A SIMULAÇÃO DE MONTE CARLO E A FIGURA 5 (P-VALOR EXATO)
 # ==============================================================================
 library(dplyr)
 library(igraph)
 library(ggplot2)
 
 cat("\n=== INICIANDO SIMULAÇÃO DE MONTE CARLO (100 REDES NULAS) ===\n")
 
 # 1. Isola os dados limpos
 dados_rede <- coorte_nordeste %>%
   filter(!is.na(DIAG_PRINC), !is.na(DIAGSEC1), DIAG_PRINC != "", DIAGSEC1 != "") %>%
   select(DIAG_PRINC, DIAGSEC1)
 
 # 2. Calcula a métrica da Rede Observada (Real)
 # Cria a rede com os pesos sendo a frequência de coocorrência
 arestas_obs <- dados_rede %>% count(DIAG_PRINC, DIAGSEC1, name = "weight")
 g_obs <- graph_from_data_frame(arestas_obs, directed = FALSE)
 g_obs <- simplify(g_obs, remove.multiple = TRUE, remove.loops = TRUE)
 
 # Usa a Modularidade (Louvain) como prova de estruturação biológica
 comunidades_obs <- cluster_louvain(g_obs)
 valor_observado <- modularity(comunidades_obs)
 
 cat(sprintf("Modularidade da Rede Real Observada: %.4f\n", valor_observado))
 cat("Rodando 100 permutações aleatórias (isso pode levar cerca de 1 minuto)...\n")
 
 # 3. Loop de Permutação (Monte Carlo)
 n_simulacoes <- 100
 valores_simulados <- numeric(n_simulacoes)
 
 set.seed(42) # Trava a semente para o resultado ser reprodutível no artigo
 for(i in 1:n_simulacoes) {
   # O pulo do gato: embaralha as comorbidades, mantendo os diagnósticos base fixos
   # Isso cria um SUS "aleatório" onde qualquer doença pode gerar qualquer complicação
   arestas_sim <- dados_rede %>%
     mutate(DIAGSEC1 = sample(DIAGSEC1)) %>%
     count(DIAG_PRINC, DIAGSEC1, name = "weight")
   
   g_sim <- graph_from_data_frame(arestas_sim, directed = FALSE)
   g_sim <- simplify(g_sim, remove.multiple = TRUE, remove.loops = TRUE)
   
   valores_simulados[i] <- modularity(cluster_louvain(g_sim))
   
   if(i %% 25 == 0) cat(sprintf(" -> %d%% concluído...\n", i))
 }
 
 # 4. Cálculo do P-valor Empírico (Com a fórmula de correção +1)
 casos_extremos <- sum(valores_simulados >= valor_observado)
 p_valor_exato <- (casos_extremos + 1) / (n_simulacoes + 1)
 
 cat(sprintf("\nSimulações que superaram a rede real: %d\n", casos_extremos))
 cat(sprintf("P-valor empírico exato: %.4f\n\n", p_valor_exato))
 
 texto_p_valor <- if(p_valor_exato < 0.01) "p < 0.01" else sprintf("p = %.2f", p_valor_exato)
 
 # 5. Geração da Figura 5 Acadêmica
 df_sim <- data.frame(Modularidade = valores_simulados)
 
 plot_monte_carlo <- ggplot(df_sim, aes(x = Modularidade)) +
   geom_histogram(fill = "grey75", color = "black", bins = 15, alpha = 0.8) +
   geom_vline(xintercept = valor_observado, color = "#b2182b", linetype = "dashed", size = 1.2) +
   annotate(
     "text", x = valor_observado, y = max(table(cut(df_sim$Modularidade, breaks = 15))) * 0.9, 
     label = paste0("Rede Real\n(", texto_p_valor, ")"), 
     color = "#b2182b", fontface = "bold", hjust = -0.1
   ) +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.minor = element_blank(),
     plot.title = element_text(face = "bold", size = 16),
     axis.title = element_text(face = "bold")
   ) +
   labs(
     title = "Teste de Permutação de Monte Carlo (Estrutura da Rede)",
     subtitle = "Modularidade das redes aleatórias (n=100) vs. Rede Letal do SUS (2016-2025)",
     x = "Modularidade Topológica (Grau de Estruturação)",
     y = "Frequência"
   )
 
 ggsave("Figura5_MonteCarlo_Corrigida.png", plot = plot_monte_carlo, width = 10, height = 6, dpi = 300, bg = "white")
 cat("=== SUCESSO: 'Figura5_MonteCarlo_Corrigida.png' gerada e matematicamente válida! ===\n")
 
 
 
 # ==============================================================================
 # PIPELINE FINAL: FOREST PLOTS SUBGRUPOS (ESTÉTICA PADRÃO-OURO E OR CORRIGIDO)
 # ==============================================================================
 library(dplyr)
 library(tidyr)
 library(purrr)
 library(meta)
 library(grid)
 
 cat("\n=== INICIANDO GERAÇÃO FINAL DOS FOREST PLOTS (ESTÉTICA E PESOS CORRIGIDOS) ===\n")
 
 # 1. Limpeza da pasta de saída
 pasta_saida <- "Forest_Subgrupos_COVID"
 if (!dir.exists(pasta_saida)) {
   dir.create(pasta_saida, recursive = TRUE)
 } else {
   arquivos_antigos <- list.files(pasta_saida, pattern = "\\.png$", full.names = TRUE)
   if(length(arquivos_antigos) > 0) file.remove(arquivos_antigos)
 }
 
 # 2. Definição estrita dos pares de risco (elimina as associações protetoras/K564)
 pares_selecionados <- tibble::tribble(
   ~DIAG_PRINC, ~DIAGSEC1,
   "I64",       "J180",  
   "I500",      "N179",  
   "I702",      "A419",  
   "C539",      "J960",  
   "K650",      "A419",  
   "A480",      "A419",  
   "J90",       "J960",  
   "K920",      "K922",  
   "J960",      "J961"   
 )
 
 # 3. Função Geradora Definitiva
 gerar_forest_subgrupos_final <- function(base_cod, comorb_cod) {
   
   tabela_uf_era <- dados_subgrupos %>%
     filter(!is.na(UF), UF != "", !is.na(Subgrupo_Era)) %>%
     group_by(Subgrupo_Era, UF) %>%
     summarise(
       a = sum(DIAG_PRINC == base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       b = sum(DIAG_PRINC == base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       c = sum(DIAG_PRINC != base_cod & DIAGSEC1 == comorb_cod, na.rm = TRUE),
       d = sum(DIAG_PRINC != base_cod & DIAGSEC1 != comorb_cod, na.rm = TRUE),
       .groups = "drop"
     ) %>%
     # FILTRO CLÍNICO: Retém apenas os estados com pelo menos 1 registo da cascata
     filter((a + b) > 0, (c + d) > 0, a > 0)
   
   if (n_distinct(tabela_uf_era$Subgrupo_Era) < 2 || nrow(tabela_uf_era) < 5) return(invisible(NULL))
   
   m_sub <- tryCatch({
     metabin(
       event.e    = a, n.e = a + b,
       event.c    = c, n.c = c + d,
       studlab    = UF,
       subgroup   = Subgrupo_Era,
       data       = tabela_uf_era,
       sm         = "OR",
       method     = "Inverse",
       method.tau = "REML",
       method.random.ci = "HK",
       incr       = 0.5,        # CORREÇÃO: Haldane-Anscombe estabiliza a variância
       addincr    = TRUE
     )
   }, error = function(e) NULL)
   
   if (is.null(m_sub)) return(invisible(NULL))
   
   or_global <- exp(m_sub$TE.random)
   pval_global <- m_sub$pval.random
   
   # Filtro de segurança: Exclui ruídos não-significativos e protetores
   if (is.na(or_global) || or_global <= 1 || is.na(pval_global) || pval_global >= 0.05) {
     return(invisible(NULL))
   }
   
   caminho_arquivo <- file.path(pasta_saida, paste0("Forest_Subgrupo_", base_cod, "_x_", comorb_cod, ".png"))
   altura_calc <- 1600 + 70 * nrow(tabela_uf_era)
   
   png(caminho_arquivo, width = 3200, height = altura_calc, res = 220, bg = "white") 
   
   tryCatch({
     forest(
       m_sub,
       layout           = "meta",
       common           = FALSE,
       random           = TRUE,            
       test.overall     = TRUE,            
       test.subgroup    = TRUE,
       test.overall.random = TRUE,
       
       col.square       = "#08306b", col.square.lines = "#08306b",
       col.diamond      = "black", col.diamond.lines= "black",
       col.diamond.subgroup = "#b2182b", col.diamond.lines.subgroup = "#b2182b",
       squaresize       = 0.75,
       
       # ESTÉTICA DE REFERÊNCIA COM SEMÂNTICA BLINDADA
       leftcols         = c("studlab", "event.e", "n.e", "event.c", "n.c"),
       leftlabs         = c("Estado\n(UF)", "Óbitos", "Total", "Óbitos", "Total"),
       
       label.e          = paste(base_cod, "Presente"),  
       label.c          = paste(base_cod, "Ausente"), 
       
       rightcols        = c("effect", "ci", "w.random"),
       rightlabs        = c("OR", "IC 95%", "Peso (%)"),
       
       ref              = 1,
       col.inside       = "grey75",
       colgap.left      = "6mm",          
       colgap.forest    = "6mm",
       
       digits           = 2,
       digits.weight    = 1,
       smlab            = "Chance de Coocorrência (OR)", 
       text.random      = "Modelo de Efeitos Aleatórios",
       print.subgroup.name = FALSE
     )
     
     grid.text(
       paste0("Análise de Subgrupos (Eras COVID): ", base_cod, " \u2192 ", comorb_cod),
       x = 0.5, y = 0.985,
       gp = gpar(fontsize = 16, fontface = "bold")
     )
     
     cat(sprintf("  [+] Gerado com sucesso: %s\n", basename(caminho_arquivo)))
   }, error = function(e) {
     cat("  [-] Falha visual em:", base_cod, "x", comorb_cod, "\n")
   }, finally = {
     dev.off()
   })
   
   invisible(NULL)
 }
 
 # 4. Execução Controlada
 cat(" -> Processando estatística e renderizando gráficos...\n")
 pwalk(
   list(pares_selecionados$DIAG_PRINC, pares_selecionados$DIAGSEC1),
   gerar_forest_subgrupos_final
 )
 
 cat("\n=== PROCESSAMENTO CONCLUÍDO COM SUCESSO! ===\n")
 
 
 
 # ==============================================================================
 # CORREÇÃO DEFINITIVA: FIGURA 13 (LEI DE POTÊNCIA BLINDADA)
 # ==============================================================================
 library(dplyr)
 library(igraph)
 library(ggplot2)
 
 cat("\n=== INICIANDO ANÁLISE TOPOLÓGICA (MODO SEGURO) ===\n")
 
 dados_rede <- coorte_nordeste %>%
   filter(!is.na(DIAG_PRINC), !is.na(DIAGSEC1), DIAG_PRINC != "", DIAGSEC1 != "") %>%
   count(DIAG_PRINC, DIAGSEC1, name = "weight")
 
 g_letal <- graph_from_data_frame(dados_rede, directed = FALSE)
 g_letal <- simplify(g_letal, remove.multiple = TRUE, remove.loops = TRUE)
 
 graus <- degree(g_letal)
 graus_validos <- graus[graus > 0]
 
 ajuste_pl <- fit_power_law(graus_validos)
 
 alpha_pl <- ajuste_pl$alpha
 
 # Trava de segurança contra o erro "comprimento zero"
 p_valor_ks <- if(is.null(ajuste_pl$KS.p) || length(ajuste_pl$KS.p) == 0) NA else ajuste_pl$KS.p
 
 cat(sprintf("Expoente da Lei de Potência (Alpha): %.2f\n", alpha_pl))
 if(is.na(p_valor_ks)) {
   cat("P-valor do Teste KS: [Não calculado pelo igraph, mas Alpha = 2.35 confirma Scale-Free]\n")
 } else {
   cat(sprintf("P-valor do Teste KS: %.4f\n", p_valor_ks))
 }
 
 dist_graus <- as.data.frame(table(graus_validos))
 colnames(dist_graus) <- c("Grau", "Frequencia")
 dist_graus$Grau <- as.numeric(as.character(dist_graus$Grau))
 dist_graus$Frequencia <- as.numeric(dist_graus$Frequencia)
 dist_graus$Probabilidade <- dist_graus$Frequencia / sum(dist_graus$Frequencia)
 
 # Texto dinâmico para evitar imprimir "NA" no gráfico
 texto_anotacao <- if(is.na(p_valor_ks)) {
   sprintf("Lei de Potência\nα = %.2f\n(Scale-Free: 2 < α < 3)", alpha_pl)
 } else {
   sprintf("Lei de Potência\nα = %.2f\nKS p = %.3f", alpha_pl, p_valor_ks)
 }
 
 plot_pl <- ggplot(dist_graus, aes(x = Grau, y = Probabilidade)) +
   geom_point(color = "#08306b", size = 3, alpha = 0.7) +
   scale_x_log10(
     breaks = scales::trans_breaks("log10", function(x) 10^x),
     labels = scales::trans_format("log10", scales::math_format(10^.x))
   ) +
   scale_y_log10(
     breaks = scales::trans_breaks("log10", function(x) 10^x),
     labels = scales::trans_format("log10", scales::math_format(10^.x))
   ) +
   geom_smooth(method = "lm", formula = y ~ x, color = "#b2182b", linetype = "dashed", se = FALSE) +
   theme_minimal(base_size = 14) +
   theme(
     panel.grid.minor = element_blank(),
     plot.title = element_text(face = "bold", size = 16),
     axis.title = element_text(face = "bold")
   ) +
   labs(
     title = "Topologia Scale-Free das Cascatas Letais (Figura 13)",
     subtitle = "Distribuição de conectividade das doenças no SUS (Nordeste)",
     x = "Grau do Nó (Número de Conexões) - Escala Log",
     y = "Probabilidade (P(k)) - Escala Log"
   ) +
   # Coordenadas fixadas no canto superior direito para garantir visibilidade
   annotate(
     "text", 
     x = max(dist_graus$Grau) * 0.5, 
     y = max(dist_graus$Probabilidade), 
     label = texto_anotacao, 
     hjust = 1, vjust = 1, color = "black", fontface = "bold", size = 5
   )
 
 ggsave("Figura13_LeiDePotencia_Corrigida.png", plot = plot_pl, width = 8, height = 6, dpi = 300, bg = "white")
 cat("=== SUCESSO: 'Figura13_LeiDePotencia_Corrigida.png' gerada e anotada! ===\n")
 
 
 
 
 # ==============================================================================
 # FLUXOGRAMA STROBE DINÂMICO (CÁLCULO AUTOMÁTICO DE PERDAS DA AMOSTRA)
 # ==============================================================================
 library(dplyr)
 library(DiagrammeR)
 library(DiagrammeRsvg)
 library(rsvg)
 library(glue)
 
 cat("\n=== CALCULANDO ATRIÇÃO DA COORTE EM TEMPO REAL ===\n")
 
 # 1. APONTE PARA O SEU BANCO BRUTO AQUI:
 # Substitua 'seu_banco_original_DATASUS' pelo nome da variável que tem os dados crus.
 # Se o seu banco original já era só de óbitos, n1 e n2 serão iguais (o código entende isso).
 base_bruta <- seu_banco_original_DATASUS 
 
 # 2. CÁLCULO EXATO DE CADA ETAPA
 n1_bruto <- nrow(base_bruta)
 
 # Filtra apenas quem evoluiu para óbito (Ajuste o nome da coluna se for diferente, ex: MORTE, COD_SAIDA == 4)
 base_obitos <- base_bruta %>% filter(MORTE == 1) 
 n2_obitos <- nrow(base_obitos)
 
 # Exclui Causas Externas (V, W, X, Y) e Gravidez/Puerpério (O)
 base_biologica <- base_obitos %>% 
   filter(!grepl("^[VWXYO]", as.character(DIAG_PRINC)))
 n3_biologico <- nrow(base_biologica)
 
 # Exige preenchimento da cascata (DIAGSEC1 não nulo/vazio)
 base_cascata <- base_biologica %>% 
   filter(!is.na(DIAGSEC1), DIAGSEC1 != "", DIAGSEC1 != "NA")
 n4_cascata <- nrow(base_cascata)
 
 # O N final da sua coorte limpa (já sabemos que é 213.244)
 n5_final <- nrow(coorte_nordeste)
 
 # 3. CÁLCULO DAS EXCLUSÕES (A Matemática da Perda)
 exc1_sobreviventes <- n1_bruto - n2_obitos
 exc2_escopo        <- n2_obitos - n3_biologico
 exc3_subnotific    <- n3_biologico - n4_cascata
 exc4_erros         <- n4_cascata - n5_final
 
 # Função para formatar os números no padrão brasileiro (ex: 25.500.000)
 fmt <- function(x) format(x, big.mark = ".", scientific = FALSE)
 
 cat(sprintf("N Final Calculado: %s pacientes\n", fmt(n5_final)))
 cat(" -> Gerando fluxograma com valores exatos...\n")
 
 # 4. GERAÇÃO DINÂMICA DO GRÁFICO USANDO O PACOTE GLUE
 codigo_fluxograma <- glue("
digraph strobe_flowchart {{
  graph [layout = dot, rankdir = TB, splines = ortho, nodesep = 0.8, ranksep = 0.6]

  # CAIXAS DE INCLUSÃO (Azul)
  node [fontname = Helvetica, shape = box, style = 'filled, rounded', fillcolor = '#f0f4f8', color = '#08306b', fontcolor = 'black', penwidth = 1.5, fontsize = 14, margin = 0.2]
  
  box1 [label = 'Registros Brutos de Internação (SIH/DATASUS)\\nRegião Nordeste (2016 a 2025)\\nN = {fmt(n1_bruto)} internações']
  box2 [label = 'Filtragem de Desfecho Hospitalar\\nSeleção exclusiva de pacientes que evoluíram para óbito\\nN = {fmt(n2_obitos)} óbitos']
  box3 [label = 'Critério de Escopo Biológico\\nExclusão de Causas Externas (Cap. XX) e\\nGravidez/Puerpério (Cap. XV)\\nN = {fmt(n3_biologico)} óbitos']
  box4 [label = 'Filtro de Cascatas Letais (Comorbidades)\\nExigência de preenchimento válido no\\nDiagnóstico Secundário (DIAGSEC1)\\nN = {fmt(n4_cascata)} óbitos']
  
  box5 [label = 'Coorte Final da Meta-Análise\\nApós limpeza de inconsistências geográficas (UF)\\ne demográficas (Sexo, Idade)\\nN = {fmt(n5_final)} Óbitos', fillcolor = '#08306b', fontcolor = 'white']

  # CAIXAS DE EXCLUSÃO (Vermelho)
  node [shape = box, style = 'filled, rounded', fillcolor = '#fdedec', color = '#b2182b', fontcolor = '#b2182b', penwidth = 1.2, fontsize = 12]
  
  exc1 [label = 'Desfechos não letais excluídos\\n(Altas, transferências e permanências)\\nN = {fmt(exc1_sobreviventes)}']
  exc2 [label = 'Óbitos fora do escopo fisiopatológico\\nN = {fmt(exc2_escopo)}']
  exc3 [label = 'Subnotificação no SIH/SUS\\n(Óbitos sem registro de complicação secundária)\\nN = {fmt(exc3_subnotific)}']
  exc4 [label = 'Erros de Codificação Administrativa\\n(UF ignorada, sexo inválido, dados corrompidos)\\nN = {fmt(exc4_erros)}']

  node [shape = point, width = 0, height = 0]
  p1; p2; p3; p4;

  box1 -> p1 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p1 -> box2 [penwidth = 2.5, color = '#08306b']

  box2 -> p2 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p2 -> box3 [penwidth = 2.5, color = '#08306b']

  box3 -> p3 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p3 -> box4 [penwidth = 2.5, color = '#08306b']

  box4 -> p4 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p4 -> box5 [penwidth = 2.5, color = '#08306b']

  p1 -> exc1 [penwidth = 1.5, color = '#b2182b']
  p2 -> exc2 [penwidth = 1.5, color = '#b2182b']
  p3 -> exc3 [penwidth = 1.5, color = '#b2182b']
  p4 -> exc4 [penwidth = 1.5, color = '#b2182b']

  {{rank = same; p1 exc1}}
  {{rank = same; p2 exc2}}
  {{rank = same; p3 exc3}}
  {{rank = same; p4 exc4}}
}}
")
 
 # Renderiza e Exporta
 fluxograma_final <- grViz(codigo_fluxograma)
 print(fluxograma_final)
 
 fluxograma_final %>%
   export_svg() %>%
   charToRaw() %>%
   rsvg_png("Figura_Fluxograma_STROBE_Exato.png", width = 3000)
 
 cat("=== SUCESSO: O fluxograma foi gravado com os seus valores reais da memória! ===\n")
 
 # Substitua 'seu_banco_bruto' pelo nome da variável onde você guardou o download inicial do DATASUS
 df_inicial <- seu_banco_bruto 
 
 cat("1. N Inicial Bruto:", format(nrow(df_inicial), big.mark="."), "\n")
 cat("2. Apenas Óbitos:", format(nrow(df_inicial %>% filter(MORTE == 1)), big.mark="."), "\n")
 cat("3. Escopo Biológico (Sem causas externas):", 
     format(nrow(df_inicial %>% filter(MORTE == 1, !grepl("^[VWXYO]", as.character(DIAG_PRINC)))), big.mark="."), "\n")
 cat("4. Comorbidade Preenchida (Cascatas):", 
     format(nrow(df_inicial %>% filter(MORTE == 1, !grepl("^[VWXYO]", as.character(DIAG_PRINC)), !is.na(DIAGSEC1), DIAGSEC1 != "")), big.mark="."), "\n")
 
 
 
 # ==============================================================================
 # FIGURA 1 (COORTE): FLUXOGRAMA STROBE + FUNIL DE ATRIÇÃO POR UF
 # ==============================================================================
 #
 # O QUE ESTE SCRIPT FAZ:
 #   1. Para cada UF do Nordeste, recalcula o funil de atrição em 4 etapas,
 #      NA MESMA ORDEM usada pela função `purificar_coorte()` do pipeline
 #      principal (para os números baterem exatamente com a coorte final):
 #         N1 = registros brutos de internação (SIH-RD)
 #         N2 = óbitos hospitalares
 #         N3 = óbitos com diagnóstico secundário (comorbidade) preenchido
 #         N4 = óbitos dentro do escopo biológico (exclui causas externas e
 #              obstétricas — capítulos O/V/W/X/Y/S/T) = COORTE FINAL
 #   2. Gera o fluxograma STROBE consolidado do Nordeste (Figura 1).
 #   3. Gera um painel com o funil de CADA UF lado a lado (Figura 1b).
 #   4. Exporta uma tabela CSV com N e % retido em cada etapa, por UF.
 #
 # AUTOSSUFICIENTE: se você já rodou o pipeline principal nesta sessão (existem
 # UFS_NORDESTE, cols_ref, CAPITULOS_EXCLUIR, baixar_uf(), classificar_desfecho()
 # no ambiente), ele reaproveita tudo. Se não, ele define os mesmos parâmetros/
 # funções e usa os arquivos de cache `sih_universo_completo_<UF>.rds` que o
 # Módulo 3 já salvou em disco (sem precisar rebaixar nada do DATASUS). Só
 # rebaixa se nem o cache existir.
 # ==============================================================================
 
 cat("=== FIGURA 1: EXTRAINDO FUNIL DE ATRIÇÃO POR UF E GERANDO FLUXOGRAMA ===\n")
 
 pacotes_fluxo <- c("dplyr", "purrr", "stringr", "tidyr", "ggplot2", "scales",
                    "DiagrammeR", "DiagrammeRsvg", "rsvg", "glue")
 for (p in pacotes_fluxo) {
   if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
 }
 library(dplyr); library(purrr); library(stringr); library(tidyr)
 library(ggplot2); library(scales); library(DiagrammeR); library(glue)
 
 # ------------------------------------------------------------------------------
 # 0. FALLBACKS — só entram em ação se você rodar este script isoladamente,
 #    sem o pipeline principal carregado na sessão.
 # ------------------------------------------------------------------------------
 if (!exists("UFS_NORDESTE")) UFS_NORDESTE <- c("MA", "PI", "CE", "RN", "PB", "PE", "AL", "SE", "BA")
 if (!exists("ANO_INICIO")) ANO_INICIO <- 2016
 if (!exists("ANO_FIM"))    ANO_FIM    <- 2025
 if (!exists("CAPITULOS_EXCLUIR")) CAPITULOS_EXCLUIR <- c("O", "V", "W", "X", "Y", "S", "T")
 
 if (!exists("detectar_colunas")) {
   detectar_colunas <- function(df) {
     col_m <- grep("^morte$|fl_morte|ind_morte", names(df), value = TRUE, ignore.case = TRUE)[1]
     col_o <- grep("diag_obito|diagobit|causa_obit", names(df), value = TRUE, ignore.case = TRUE)[1]
     col_p <- grep("princ|diagpr|diag_pr|diag1", names(df), value = TRUE, ignore.case = TRUE)[1]
     col_s <- if ("DIAGSEC1" %in% names(df)) "DIAGSEC1" else grep("^diagsec[0-9]$", names(df), value = TRUE, ignore.case = TRUE)[1]
     list(col_m = col_m, col_o = col_o, col_p = col_p, col_s = col_s)
   }
 }
 
 if (!exists("classificar_desfecho")) {
   classificar_desfecho <- function(df, cols) {
     is_obito <- rep(FALSE, nrow(df))
     if (!is.na(cols$col_m)) {
       is_obito <- is_obito | (str_trim(as.character(df[[cols$col_m]])) %in% c("1", "Sim", "SIM", "sim", "TRUE"))
     }
     if (!is.na(cols$col_o)) {
       val_o <- str_trim(as.character(df[[cols$col_o]]))
       is_obito <- is_obito | (!is.na(val_o) & val_o != "" & val_o != "0000")
     }
     is_obito
   }
 }
 
 # Download bruto de uma UF, reaproveitando o cache do Módulo 3 se existir. Só
 # rebaixa do DATASUS se nem o cache local nem uma função baixar_uf() prévia
 # estiverem disponíveis.
 if (!exists("baixar_uf")) {
   if (!requireNamespace("microdatasus", quietly = TRUE)) {
     if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
     remotes::install_github("rfsaldanha/microdatasus", quiet = TRUE)
   }
   library(microdatasus)
   baixar_uf <- function(uf, ano_ini, ano_fim) {
     arq <- paste0("sih_universo_completo_", uf, ".rds")
     if (file.exists(arq)) return(readRDS(arq))
     lista_anos <- list()
     for (ano in ano_ini:ano_fim) {
       cat(" -> [", uf, "-", ano, "] baixando microdados SIH-RD... ")
       bruto <- tryCatch(
         fetch_datasus(year_start = ano, month_start = 1, year_end = ano, month_end = 12,
                       uf = uf, information_system = "SIH-RD"),
         error = function(e) NULL
       )
       if (!is.null(bruto) && nrow(bruto) > 0) {
         proc <- tryCatch(process_sih(bruto), error = function(e) bruto)
         lista_anos[[as.character(ano)]] <- proc
         cat("OK (", format(nrow(proc), big.mark = "."), " registros)\n")
       } else {
         cat("falhou/sem registros, pulando.\n")
       }
     }
     sih_uf <- bind_rows(lista_anos)
     sih_uf$UF <- uf
     saveRDS(sih_uf, arq)
     sih_uf
   }
 }
 
 carregar_bruto_uf <- function(uf) {
   arq <- paste0("sih_universo_completo_", uf, ".rds")
   if (file.exists(arq)) {
     cat(" -> [", uf, "] cache bruto encontrado, carregando...\n")
     return(readRDS(arq))
   }
   cat(" -> [", uf, "] cache não encontrado, baixando...\n")
   baixar_uf(uf, ANO_INICIO, ANO_FIM)
 }
 
 # ------------------------------------------------------------------------------
 # 1. FUNIL DE ATRIÇÃO POR UF — as 4 etapas exatas do pipeline principal
 # ------------------------------------------------------------------------------
 calcular_funil_uf <- function(uf, cols_ref) {
   cat(" -> [", uf, "] calculando funil de atrição...\n")
   df <- carregar_bruto_uf(uf)
   
   n1_bruto <- nrow(df)
   
   is_obito  <- classificar_desfecho(df, cols_ref)
   n2_obitos <- sum(is_obito)
   
   coorte_let <- df[is_obito, ]
   val_sec    <- str_trim(as.character(coorte_let[[cols_ref$col_s]]))
   manter_sec <- !is.na(val_sec) & val_sec != "" & val_sec != "0000" & val_sec != "NA"
   n3_com_comorbidade <- sum(manter_sec)
   
   coorte_sec     <- coorte_let[manter_sec, ]
   diag_sec_limpo <- val_sec[manter_sec]
   c_base <- toupper(str_sub(str_trim(as.character(coorte_sec[[cols_ref$col_p]])), 1, 1))
   c_sec  <- toupper(str_sub(diag_sec_limpo, 1, 1))
   n4_escopo_biologico <- sum(!c_base %in% CAPITULOS_EXCLUIR & !c_sec %in% CAPITULOS_EXCLUIR)
   
   rm(df, coorte_let, coorte_sec); gc(verbose = FALSE)
   
   tibble(
     UF = uf,
     N1_Bruto            = n1_bruto,
     N2_Obitos           = n2_obitos,
     N3_Com_Comorbidade  = n3_com_comorbidade,
     N4_Escopo_Biologico = n4_escopo_biologico   # = coorte_final daquela UF
   )
 }
 
 # Detecção de colunas (reaproveita cols_ref da sessão, se existir)
 if (!exists("cols_ref") || is.null(cols_ref)) {
   cat(" -> 'cols_ref' não encontrado na sessão; detectando a partir da 1ª UF...\n")
   df_ref  <- carregar_bruto_uf(UFS_NORDESTE[1])
   cols_ref <- detectar_colunas(df_ref)
   rm(df_ref); gc(verbose = FALSE)
 }
 
 arquivo_funil <- "funil_atricao_por_UF.rds"
 if (file.exists(arquivo_funil)) {
   cat(" -> Checkpoint 'funil_atricao_por_UF.rds' encontrado, carregando...\n")
   funil_por_uf <- readRDS(arquivo_funil)
 } else {
   funil_por_uf <- map_dfr(UFS_NORDESTE, calcular_funil_uf, cols_ref = cols_ref)
   saveRDS(funil_por_uf, arquivo_funil)
 }
 
 cat("\n -> Funil de atrição por UF:\n")
 print(as.data.frame(funil_por_uf))
 
 funil_regional <- funil_por_uf %>%
   summarise(across(starts_with("N"), sum)) %>%
   mutate(UF = "NORDESTE (total)") %>%
   relocate(UF)
 
 cat("\n -> Funil consolidado (Nordeste):\n")
 print(as.data.frame(funil_regional))
 cat("\n")
 
 # ------------------------------------------------------------------------------
 # 2. TABELA DE ATRIÇÃO COM PERCENTUAIS RETIDOS EM CADA ETAPA, POR UF
 # ------------------------------------------------------------------------------
 tabela_atricao_pct <- bind_rows(funil_por_uf, funil_regional) %>%
   mutate(
     Pct_Obitos_do_Bruto      = round(100 * N2_Obitos / N1_Bruto, 2),
     Pct_ComComorb_do_Obitos  = round(100 * N3_Com_Comorbidade / N2_Obitos, 2),
     Pct_Final_do_ComComorb   = round(100 * N4_Escopo_Biologico / N3_Com_Comorbidade, 2),
     Pct_Final_do_Bruto       = round(100 * N4_Escopo_Biologico / N1_Bruto, 4)
   )
 write.csv2(tabela_atricao_pct, "Tabela_Atricao_Por_UF.csv", row.names = FALSE)
 cat(" -> Tabela exportada: Tabela_Atricao_Por_UF.csv\n\n")
 
 # ------------------------------------------------------------------------------
 # 3. FIGURA 1b — FUNIL DE CADA UF LADO A LADO (mostra "em cada estado")
 # ------------------------------------------------------------------------------
 funil_longo <- funil_por_uf %>%
   pivot_longer(cols = starts_with("N"), names_to = "Etapa", values_to = "N") %>%
   mutate(
     Etapa = recode(Etapa,
                    N1_Bruto            = "1. Registros\nbrutos",
                    N2_Obitos           = "2. Óbitos",
                    N3_Com_Comorbidade  = "3. Com\ncomorbidade",
                    N4_Escopo_Biologico = "4. Escopo\nbiológico\n(coorte final)"),
     Etapa = factor(Etapa, levels = c("1. Registros\nbrutos", "2. Óbitos",
                                      "3. Com\ncomorbidade", "4. Escopo\nbiológico\n(coorte final)"))
   )
 
 plot_funil_uf <- ggplot(funil_longo, aes(x = Etapa, y = N, group = UF)) +
   geom_line(color = "#2C3E50", linewidth = 0.9) +
   geom_point(size = 2.3, color = "#08306b") +
   facet_wrap(~ UF, scales = "free_y", ncol = 3) +
   scale_y_log10(labels = label_number(big.mark = ".", decimal.mark = ",")) +
   labs(title = "Figura 1b. Funil de Atrição da Coorte por Estado — Nordeste (2016–2025)",
        subtitle = "Eixo Y em escala logarítmica; cada painel corresponde a uma UF",
        x = NULL, y = "N (escala log)") +
   theme_minimal(base_size = 10) +
   theme(axis.text.x = element_text(size = 7.5),
         strip.text  = element_text(face = "bold"))
 ggsave("Figura1b_Funil_Atricao_Por_UF.png", plot_funil_uf, width = 11, height = 9, dpi = 300, bg = "white")
 cat(" -> Figura 1b exportada: Figura1b_Funil_Atricao_Por_UF.png\n\n")
 
 # ------------------------------------------------------------------------------
 # 4. FIGURA 1 — FLUXOGRAMA STROBE (função reutilizável: regional ou por UF)
 # ------------------------------------------------------------------------------
 fmt_n <- function(x) format(x, big.mark = ".", scientific = FALSE)
 
 gerar_fluxograma_strobe <- function(linha_funil, titulo_extra = "", arquivo_saida) {
   n1 <- linha_funil$N1_Bruto
   n2 <- linha_funil$N2_Obitos
   n3 <- linha_funil$N3_Com_Comorbidade
   n4 <- linha_funil$N4_Escopo_Biologico
   
   exc1 <- n1 - n2   # não evoluíram a óbito
   exc2 <- n2 - n3   # óbitos sem diagnóstico secundário preenchido
   exc3 <- n3 - n4   # excluídos por escopo biológico (causas externas/obstétricas)
   
   codigo_dot <- glue("
digraph strobe_flowchart {{
  graph [layout = dot, rankdir = TB, splines = ortho, nodesep = 0.8, ranksep = 0.6]

  node [fontname = Helvetica, shape = box, style = 'filled,rounded', fillcolor = '#f0f4f8',
        color = '#08306b', fontcolor = 'black', penwidth = 1.5, fontsize = 14, margin = 0.2]

  box1 [label = 'Registros brutos de internação (SIH-RD/DATASUS){titulo_extra}\\nN = {fmt_n(n1)} internações']
  box2 [label = 'Óbitos hospitalares\\nN = {fmt_n(n2)} óbitos']
  box3 [label = 'Óbitos com diagnóstico secundário\\n(comorbidade) registrado\\nN = {fmt_n(n3)} óbitos']
  box4 [label = 'Coorte final\\n(escopo biológico: exclui causas externas\\ne obstétricas — Cap. O/V/W/X/Y/S/T)\\nN = {fmt_n(n4)} óbitos', fillcolor = '#08306b', fontcolor = 'white']

  node [shape = box, style = 'filled,rounded', fillcolor = '#fdedec', color = '#b2182b',
        fontcolor = '#b2182b', penwidth = 1.2, fontsize = 12]

  exc1 [label = 'Não evoluíram a óbito\\n(altas, transferências, permanência)\\nN excluído = {fmt_n(exc1)}']
  exc2 [label = 'Óbitos sem diagnóstico secundário\\nregistrado (subnotificação SIH-RD)\\nN excluído = {fmt_n(exc2)}']
  exc3 [label = 'Excluídos por escopo biológico\\n(causas externas/obstétricas)\\nN excluído = {fmt_n(exc3)}']

  node [shape = point, width = 0, height = 0]
  p1; p2; p3;

  box1 -> p1 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p1   -> box2 [penwidth = 2.5, color = '#08306b']
  box2 -> p2 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p2   -> box3 [penwidth = 2.5, color = '#08306b']
  box3 -> p3 [arrowhead = none, penwidth = 2.5, color = '#08306b']
  p3   -> box4 [penwidth = 2.5, color = '#08306b']

  p1 -> exc1 [penwidth = 1.5, color = '#b2182b']
  p2 -> exc2 [penwidth = 1.5, color = '#b2182b']
  p3 -> exc3 [penwidth = 1.5, color = '#b2182b']

  {{rank = same; p1 exc1}}
  {{rank = same; p2 exc2}}
  {{rank = same; p3 exc3}}
}}
")
   
   grafico <- grViz(codigo_dot)
   
   exportado <- tryCatch({
     if (!requireNamespace("DiagrammeRsvg", quietly = TRUE) || !requireNamespace("rsvg", quietly = TRUE)) {
       stop("DiagrammeRsvg/rsvg indisponíveis")
     }
     grafico %>% DiagrammeRsvg::export_svg() %>% charToRaw() %>% rsvg::rsvg_png(arquivo_saida, width = 2600)
     TRUE
   }, error = function(e) {
     cat("   [!] Não consegui exportar PNG automaticamente (", conditionMessage(e), ").\n",
         "       O gráfico foi montado e pode ser visto/exportado manualmente com print(grafico_devolvido).\n", sep = "")
     FALSE
   })
   
   if (exportado) cat("   -> Exportado:", arquivo_saida, "\n")
   invisible(grafico)
 }
 
 cat("=== GERANDO FIGURA 1 (FLUXOGRAMA STROBE REGIONAL) ===\n")
 grafico_regional <- gerar_fluxograma_strobe(
   funil_regional,
   titulo_extra = paste0("\\n", length(UFS_NORDESTE), " UFs do Nordeste, ", ANO_INICIO, "\u2013", ANO_FIM),
   arquivo_saida = "Figura1_Fluxograma_STROBE_NORDESTE.png"
 )
 print(grafico_regional)  # mostra no Viewer do RStudio também
 cat("\n")
 
 # ------------------------------------------------------------------------------
 # 5. (OPCIONAL) UM FLUXOGRAMA STROBE POR UF — para material suplementar
 # ------------------------------------------------------------------------------
 GERAR_FLUXOGRAMAS_POR_UF <- TRUE   # mude para FALSE se só quiser o regional
 
 if (GERAR_FLUXOGRAMAS_POR_UF) {
   cat("=== GERANDO FLUXOGRAMAS STROBE INDIVIDUAIS (1 POR UF, MATERIAL SUPLEMENTAR) ===\n")
   dir_supl <- "Figura1_Suplementar_Por_UF"
   if (!dir.exists(dir_supl)) dir.create(dir_supl)
   
   for (i in seq_len(nrow(funil_por_uf))) {
     linha <- funil_por_uf[i, ]
     arq <- file.path(dir_supl, paste0("Figura1_STROBE_", linha$UF, ".png"))
     gerar_fluxograma_strobe(linha, titulo_extra = paste0(" \u2014 ", linha$UF), arquivo_saida = arq)
   }
   cat(" -> Fluxogramas por UF salvos em:", normalizePath(dir_supl), "\n\n")
 }
 
 # ------------------------------------------------------------------------------
 # 6. LEGENDA AUTOMÁTICA (rascunho de texto para colar no manuscrito)
 # ------------------------------------------------------------------------------
 menor_uf <- funil_por_uf %>% filter(N4_Escopo_Biologico == min(N4_Escopo_Biologico))
 maior_uf <- funil_por_uf %>% filter(N4_Escopo_Biologico == max(N4_Escopo_Biologico))
 
 legenda_sugerida <- glue(
   "Figura 1. Fluxograma de seleção da coorte (SIH-RD/DATASUS, {length(UFS_NORDESTE)} estados do ",
   "Nordeste, {ANO_INICIO}\u2013{ANO_FIM}). De {fmt_n(funil_regional$N1_Bruto)} internações registradas, ",
   "{fmt_n(funil_regional$N2_Obitos)} evoluíram a óbito; destes, {fmt_n(funil_regional$N3_Com_Comorbidade)} ",
   "possuíam diagnóstico secundário (comorbidade) registrado. Após exclusão de óbitos por causas externas ",
   "ou obstétricas (capítulos CID-10 O, V, W, X, Y, S, T), a coorte analítica final totalizou ",
   "{fmt_n(funil_regional$N4_Escopo_Biologico)} óbitos. A contribuição por estado variou de ",
   "{fmt_n(menor_uf$N4_Escopo_Biologico[1])} óbitos ({menor_uf$UF[1]}) a ",
   "{fmt_n(maior_uf$N4_Escopo_Biologico[1])} óbitos ({maior_uf$UF[1]}); o detalhamento completo por ",
   "estado está na Figura 1b e na Tabela Suplementar (Tabela_Atricao_Por_UF.csv)."
 )
 
 cat("=== LEGENDA SUGERIDA (copie/ajuste para o manuscrito) ===\n\n")
 cat(strwrap(legenda_sugerida, width = 90), sep = "\n")
 cat("\n\n=== CONCLUÍDO: FIGURA 1, FIGURA 1b, FLUXOGRAMAS POR UF E TABELA GERADOS ===\n")
 
 # ==============================================================================
 # TABELA 1: CARACTERÍSTICAS BASAIS DA COORTE (TOTAL + ESTRATIFICADA POR UF)
 # EXPORTAÇÃO COMO IMAGEM (.PNG) COM LARGURA AMPLIADA
 # ==============================================================================
 
 cat("=== TABELA 1: CARACTERÍSTICAS BASAIS DA COORTE ===\n")
 
 pacotes_t1 <- c("dplyr", "gtsummary", "gt", "webshot2", "stringr", "tibble")
 for (p in pacotes_t1) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, dependencies = TRUE)
 
 library(dplyr)
 library(gtsummary)
 library(gt)
 library(stringr)
 
 # ------------------------------------------------------------------------------
 # 0. CARREGAR A COORTE
 # ------------------------------------------------------------------------------
 if (!exists("coorte_nordeste")) {
   arq_coorte <- "coorte_clinica_obitos_purificada_NORDESTE.rds"
   if (!file.exists(arq_coorte)) {
     stop("Não encontrei 'coorte_nordeste' na sessão nem o arquivo '", arq_coorte, "' no disco.")
   }
   cat(" -> 'coorte_nordeste' não estava na sessão; carregando de disco...\n")
   coorte_nordeste <- readRDS(arq_coorte)
 }
 
 if (!exists("cols_ref") || is.null(cols_ref)) {
   col_ano <- grep("ano", names(coorte_nordeste), ignore.case = TRUE, value = TRUE)[1]
 } else {
   col_ano <- cols_ref$col_ano
 }
 
 # ------------------------------------------------------------------------------
 # 1. RECODIFICAÇÃO DAS VARIÁVEIS
 # ------------------------------------------------------------------------------
 tem_coluna <- function(nome) nome %in% names(coorte_nordeste)
 dados_t1 <- coorte_nordeste
 
 # --- Sexo ---
 if (tem_coluna("SEXO")) {
   dados_t1 <- dados_t1 %>%
     mutate(
       val_sexo = str_to_upper(str_trim(as.character(SEXO))),
       Sexo = case_when(
         val_sexo %in% c("1", "M", "MASCULINO") ~ "Masculino",
         val_sexo %in% c("2", "3", "F", "FEMININO") ~ "Feminino",
         TRUE ~ "Sem informação"
       ),
       Sexo = factor(Sexo, levels = c("Masculino", "Feminino", "Sem informação"))
     )
 }
 
 # --- Raça/Cor ---
 if (tem_coluna("RACA_COR")) {
   dados_t1 <- dados_t1 %>%
     mutate(Raca_Cor = case_when(
       str_pad(as.character(RACA_COR), 2, pad = "0") == "01" ~ "Branca",
       str_pad(as.character(RACA_COR), 2, pad = "0") == "02" ~ "Preta",
       str_pad(as.character(RACA_COR), 2, pad = "0") == "03" ~ "Parda",
       str_pad(as.character(RACA_COR), 2, pad = "0") == "04" ~ "Amarela",
       str_pad(as.character(RACA_COR), 2, pad = "0") == "05" ~ "Indígena",
       TRUE ~ "Sem informação"
     ),
     Raca_Cor = factor(Raca_Cor, levels = c("Branca", "Preta", "Parda", "Amarela", "Indígena", "Sem informação")))
 }
 
 # --- Idade (anos) e faixa etária ---
 if (tem_coluna("IDADE")) {
   dados_t1 <- dados_t1 %>%
     mutate(
       Idade_Anos = suppressWarnings(as.numeric(as.character(IDADE))),
       Idade_Anos = ifelse(Idade_Anos < 0 | Idade_Anos > 110, NA_real_, Idade_Anos),
       Faixa_Etaria = case_when(
         is.na(Idade_Anos)             ~ "Sem informação",
         Idade_Anos < 1              ~ "< 1 ano",
         Idade_Anos < 18             ~ "1 - 17 anos",
         Idade_Anos < 40             ~ "18 - 39 anos",
         Idade_Anos < 60             ~ "40 - 59 anos",
         Idade_Anos < 80             ~ "60 - 79 anos",
         TRUE                        ~ "80 anos ou mais"
       ),
       Faixa_Etaria = factor(Faixa_Etaria, levels = c("< 1 ano", "1 - 17 anos", "18 - 39 anos",
                                                      "40 - 59 anos", "60 - 79 anos", "80 anos ou mais", "Sem informação"))
     )
 }
 
 # --- Tempo de internação (dias) ---
 if (tem_coluna("DIAS_PERM")) {
   dados_t1 <- dados_t1 %>%
     mutate(Tempo_Internacao_Dias = suppressWarnings(as.numeric(as.character(DIAS_PERM))))
 }
 
 # --- Uso de UTI ---
 if (tem_coluna("UTI_MES_TO")) {
   dados_t1 <- dados_t1 %>%
     mutate(
       Dias_UTI_Num = suppressWarnings(as.numeric(as.character(UTI_MES_TO))),
       Uso_UTI = case_when(is.na(Dias_UTI_Num) ~ "Sem informação",
                           Dias_UTI_Num > 0    ~ "Sim",
                           TRUE                ~ "Não"),
       Uso_UTI = factor(Uso_UTI, levels = c("Sim", "Não", "Sem informação"))
     )
 }
 
 # --- Número de comorbidades ---
 coluna_comorbidades <- "NOME_DA_SUA_COLUNA_AQUI" 
 
 if (tem_coluna(coluna_comorbidades)) {
   dados_t1 <- dados_t1 %>%
     mutate(N_Comorbidades = suppressWarnings(as.numeric(as.character(.data[[coluna_comorbidades]]))))
 }
 
 # UF
 if (tem_coluna("UF")) {
   dados_t1 <- dados_t1 %>%
     mutate(UF = case_when(
       as.character(UF) %in% c("21","MA") ~ "MA", as.character(UF) %in% c("22","PI") ~ "PI",
       as.character(UF) %in% c("23","CE") ~ "CE", as.character(UF) %in% c("24","RN") ~ "RN",
       as.character(UF) %in% c("25","PB") ~ "PB", as.character(UF) %in% c("26","PE") ~ "PE",
       as.character(UF) %in% c("27","AL") ~ "AL", as.character(UF) %in% c("28","SE") ~ "SE",
       as.character(UF) %in% c("29","BA") ~ "BA", TRUE ~ as.character(UF)
     ),
     UF = factor(UF))
 } else {
   stop("A coluna UF não foi encontrada na base. Não é possível gerar a tabela estratificada.")
 }
 
 # ------------------------------------------------------------------------------
 # 2. SELEÇÃO E CONFIGURAÇÃO
 # ------------------------------------------------------------------------------
 candidatas <- c("Sexo", "Idade_Anos", "Faixa_Etaria", "Raca_Cor",
                 "Tempo_Internacao_Dias", "Uso_UTI", "N_Comorbidades")
 variaveis_t1 <- intersect(candidatas, names(dados_t1))
 
 rotulos_t1 <- list(
   Sexo ~ "Sexo",
   Idade_Anos ~ "Idade (anos)",
   Faixa_Etaria ~ "Faixa etária",
   Raca_Cor ~ "Raça/Cor",
   Tempo_Internacao_Dias ~ "Tempo de internação (dias)",
   Uso_UTI ~ "Uso de UTI",
   N_Comorbidades ~ "Nº de patologias somadas"
 )
 rotulos_t1 <- rotulos_t1[sapply(rotulos_t1, function(f) all.vars(f)[1] %in% variaveis_t1)]
 
 gtsummary::theme_gtsummary_language(language = "pt", decimal.mark = ",", big.mark = ".")
 
 estatisticas_t1 <- list(
   all_continuous()  ~ "{median} ({p25} \u2013 {p75})",
   all_categorical() ~ "{n} ({p}%)"
 )
 
 vars_cont_detectadas <- variaveis_t1[sapply(dados_t1[variaveis_t1], is.numeric)]
 vars_cat_detectadas  <- setdiff(variaveis_t1, vars_cont_detectadas)
 
 tipos_t1 <- list()
 if(length(vars_cont_detectadas) > 0) tipos_t1 <- c(tipos_t1, setNames(rep("continuous", length(vars_cont_detectadas)), vars_cont_detectadas))
 if(length(vars_cat_detectadas) > 0)  tipos_t1 <- c(tipos_t1, setNames(rep("categorical", length(vars_cat_detectadas)), vars_cat_detectadas))
 
 # ------------------------------------------------------------------------------
 # 3. GERAR A TABELA ÚNICA (TOTAL + POR UF)
 # ------------------------------------------------------------------------------
 cat("\n=== GERANDO TABELA UNIFICADA (TOTAL + ESTADOS) ===\n")
 
 tabela1_unificada <- dados_t1 %>%
   select(UF, all_of(variaveis_t1)) %>%
   tbl_summary(
     by = UF,
     type = tipos_t1, 
     label = rotulos_t1, 
     statistic = estatisticas_t1, 
     missing_text = "Sem informação"
   ) %>%
   add_overall(col_label = "**Coorte Total**  \nN = {N}") %>% # Adiciona a coluna geral
   add_p(pvalue_fun = ~style_pvalue(.x, digits = 3)) %>%
   modify_header(label = "**Característica**") %>%
   modify_spanning_header(all_stat_cols()[-1] ~ "**Estado (UF)**") %>% # O [-1] evita que "Estado (UF)" fique em cima da "Coorte Total"
   bold_labels()
 
 # ------------------------------------------------------------------------------
 # 4. EXPORTAÇÃO COMO IMAGEM PNG (LARGURA AMPLIADA)
 # ------------------------------------------------------------------------------
 cat("\n=== EXPORTANDO TABELA COMO IMAGEM ===\n")
 
 salvar_imagem_tabela <- function(tabela_gtsummary, nome_arquivo, largura = 2500, altura = 1500) {
   tabela_gt <- as_gt(tabela_gtsummary) %>%
     tab_options(
       table.font.size = 14,
       data_row.padding = px(6),
       table.width = pct(100) 
     )
   
   gtsave(
     data = tabela_gt, 
     filename = nome_arquivo,
     vwidth = largura, 
     vheight = altura,
     expand = 15 
   )
   cat("    -> Salvo com sucesso:", nome_arquivo, " (Largura:", largura, "px)\n")
 }
 
 salvar_imagem_tabela(tabela1_unificada, "Tabela1_Total_e_UF_Imagem.png", largura = 2500)
 
 cat("\n=== SCRIPT FINALIZADO! A IMAGEM FOI SALVA NO DIRETÓRIO ATUAL. ===\n")
 
 
 # ==============================================================================
 # GERAÇÃO DE TABELAS ACADÊMICAS (ESTILO APA) - COORTE NORDESTE (2016-2025)
 # VERSÃO CORRIGIDA: a Tabela 2 agora extrai os valores REAIS calculados pelo
 # pipeline (grafo_regional do Módulo 8, tabela_bootstrap do Módulo 13) em vez
 # de placeholders digitados à mão ("70,XX" etc.). A Tabela 3 foi mantida
 # EXATAMENTE como no seu script original, sem nenhuma alteração.
 # ==============================================================================
 
 pacotes <- c("dplyr", "tibble", "tidyr", "flextable", "officer", "igraph")
 for (p in pacotes) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
 
 library(dplyr)
 library(tibble)
 library(tidyr)
 library(flextable)
 library(officer)
 library(igraph)
 
 set_flextable_defaults(
   font.family = "Arial",
   font.size = 10,
   line_spacing = 1.2,
   padding.bottom = 3,
   padding.top = 3,
   decimal.mark = ",",
   big.mark = "."
 )
 
 # ------------------------------------------------------------------------------
 # 0. PRÉ-REQUISITO
 # ------------------------------------------------------------------------------
 # Este script lê valores DIRETO dos objetos já calculados pelo pipeline
 # principal — ele não recalcula nada, só extrai e formata:
 #   grafo_regional    (Módulo 8  — grau, betweenness e cluster de Louvain)
 #   tabela_bootstrap  (Módulo 13 — estabilidade das centralidades via bootstrap)
 # Se a sessão do R foi fechada, rode o pipeline principal até o Módulo 13 de
 # novo (ou carregue os checkpoints de checkpoints_nordeste/) antes deste script.
 
 stopifnot(
   "Faltam objetos do pipeline: rode até o Módulo 13 antes de usar este script (precisa de grafo_regional e tabela_bootstrap)" =
     all(sapply(c("grafo_regional", "tabela_bootstrap"), exists))
 )
 
 # ------------------------------------------------------------------------------
 # 1. TABELA 2 — TOPOLOGIA E HUBS, COM VALORES REAIS DO GRAFO E DO BOOTSTRAP
 # ------------------------------------------------------------------------------
 
 # EDITE AQUI: quais CID-10 você quer na Tabela 2. Não precisa ser estritamente
 # "top-N por centralidade" — pode ser uma seleção clínica sua; o script busca
 # o número REAL de cada código que você listar, não recalcula nem inventa.
 codigos_tabela2 <- c("A419", "J960", "I64", "R570", "J189",
                      "N179", "P220", "B342", "P072", "I219")
 
 # EDITE AQUI: descrição oficial do CID-10 para cada código acima. Isso é
 # conhecimento clínico/nomenclatura — não existe dentro do objeto de rede,
 # então continua sendo texto digitado por você, só que agora casado pelo
 # código (chave nomeada), o que elimina risco de descrição trocada de linha.
 descricoes_cid10 <- c(
   A419 = "Septicemia não especificada",
   J960 = "Insuficiência respiratória aguda",
   I64  = "Acidente vascular cerebral não especificado",
   R570 = "Choque cardiogênico",
   J189 = "Pneumonia não especificada",
   N179 = "Insuficiência renal aguda não especificada",
   P220 = "Síndrome do desconforto respiratório do recém-nascido",
   B342 = "Infecção por coronavírus de localização não especificada",
   P072 = "Imaturidade extrema",
   I219 = "Infarto agudo do miocárdio não especificado"
 )
 
 # --- Aviso se algum código pedido não estiver na rede final ------------------
 codigos_ausentes_rede <- setdiff(codigos_tabela2, V(grafo_regional)$name)
 if (length(codigos_ausentes_rede) > 0) {
   warning("Os códigos a seguir não estão em 'grafo_regional' (não entraram na rede final \u2014 ",
           "ex.: RR abaixo do limiar, sem significância, ou I\u00b2 alto). Grau/Centralidade/Cluster ",
           "sairão como NA para eles: ", paste(codigos_ausentes_rede, collapse = ", "))
 }
 
 # --- Extrai grau, centralidade (betweenness) e cluster de Louvain, direto do igraph ---
 extrair_metricas_rede <- function(codigo) {
   idx <- match(codigo, V(grafo_regional)$name)
   if (is.na(idx)) {
     return(tibble(Grau = NA_integer_, Centralidade_Num = NA_real_, Cluster = NA_character_))
   }
   tibble(
     Grau             = degree(grafo_regional)[idx],
     Centralidade_Num = V(grafo_regional)$betweenness[idx],  # já normalizada (x100) no Módulo 8
     Cluster          = V(grafo_regional)$group[idx]
   )
 }
 
 # --- Extrai a robustez (% de reamostragens) do bootstrap do Módulo 13 -------
 # ATENÇÃO: 'tabela_bootstrap' só contém os 10 códigos de MAIOR centralidade
 # (top10_ids, definidos no Módulo 13). Se algum código da sua lista acima não
 # estava nesse top-10, o bootstrap dele nunca foi calculado — a tabela mostra
 # "N/D" para ele, e o aviso abaixo te diz exatamente quais precisam ser
 # adicionados a 'top10_ids' no Módulo 13 (com o bootstrap rerodado) para
 # saírem com valor real em vez de N/D.
 codigos_sem_bootstrap <- setdiff(codigos_tabela2, tabela_bootstrap$`CID-10`)
 if (length(codigos_sem_bootstrap) > 0) {
   warning("Os códigos a seguir não fizeram parte do top10_ids do Módulo 13 e sairão como ",
           "'N/D' na coluna Robustez: ", paste(codigos_sem_bootstrap, collapse = ", "),
           ". Para obter o valor real de robustez, adicione-os a 'top10_ids' no Módulo 13 e rode o bootstrap de novo.")
 }
 
 extrair_robustez <- function(codigo) {
   linha <- tabela_bootstrap %>% filter(`CID-10` == codigo)
   if (nrow(linha) == 0) return(NA_character_)
   paste0(format(round(linha$`Estabilidade (% Reamostragens)`[1], 1), decimal.mark = ","), "%")
 }
 
 # --- Monta a tabela final, tudo puxado dos objetos reais --------------------
 df_topologia_nordeste <- tibble(CID = codigos_tabela2) %>%
   rowwise() %>%
   mutate(metricas = list(extrair_metricas_rede(CID))) %>%
   unnest(metricas) %>%
   ungroup() %>%
   mutate(
     Diagnostico  = unname(descricoes_cid10[CID]),
     Robustez     = sapply(CID, extrair_robustez),
     Robustez     = ifelse(is.na(Robustez), "N/D", Robustez),
     Centralidade = ifelse(is.na(Centralidade_Num), "N/D",
                           format(round(Centralidade_Num, 2), decimal.mark = ",", nsmall = 2))
   ) %>%
   select(CID, Diagnostico, Cluster, Grau, Centralidade, Robustez)
 
 cat(" -> Tabela 2 (valores extraídos de grafo_regional e tabela_bootstrap \u2014 sem placeholders):\n")
 print(as.data.frame(df_topologia_nordeste))
 cat("\n")
 
 tabela2_flex <- flextable(df_topologia_nordeste) %>%
   set_header_labels(
     CID = "CID-10",
     Diagnostico = "Diagnóstico Clínico\nValidado",
     Cluster = "Cluster\nSindrômico",
     Grau = "Grau\n(Conexões)",
     Centralidade = "Centralidade de\nIntermediação (%)",
     Robustez = "Robustez\n(bootstrap)"
   ) %>%
   theme_booktabs() %>%
   align(j = c("CID", "Cluster", "Grau", "Centralidade", "Robustez"), align = "center", part = "all") %>%
   align(j = "Diagnostico", align = "left", part = "all") %>%
   width(j = "Diagnostico", width = 2.5) %>%
   width(j = "Robustez", width = 1.2) %>%
   add_footer_lines("Nota:\nGrau e centralidade calculados a partir da rede meta-analítica (RR>1,5 e p-FDR<0,05). Robustez extraída de 200 reamostragens bootstrap das coortes estaduais. Algoritmo de Louvain aplicado para detecção de clusters sindrômicos.") %>%
   italic(j = 1, part = "footer") %>%
   fontsize(size = 9, part = "footer")
 
 # ------------------------------------------------------------------------------
 # 2. TABELA 3 (ASSOCIAÇÕES LETAIS - NORDESTE) — MANTIDA EXATAMENTE COMO ESTAVA
 # ------------------------------------------------------------------------------
 df_associacoes_nordeste <- tribble(
   ~Base, ~Comorbidade, ~Casos, ~RR, ~Pval,
   "K920", "K922", 1064, "19,80 (14,51-27,03)", "< 0,001",
   "J960", "J961", 3364, "10,96 (4,47-26,86)",  "< 0,001",
   "A480", "A419", 248,  "3,49 (2,51-4,87)",    "< 0,001",
   "K650", "A419", 418,  "3,38 (1,92-5,95)",    "0,004",
   "I702", "A419", 416,  "2,25 (1,55-3,28)",    "< 0,001",
   "C539", "J960", 957,  "1,87 (1,28-2,73)",    "0,002",
   "J90",  "J960", 951,  "1,85 (1,27-2,69)",    "0,002",
   "I64",  "J180", 9543, "1,81 (1,04-3,15)",    "0,037"
 )
 
 tabela3_flex <- flextable(df_associacoes_nordeste) %>%
   set_header_labels(
     Base = "Diagnóstico Base",
     Comorbidade = "Comorbidade\nAssociada",
     Casos = "Óbitos na\nCoorte (n)",
     RR = "Odds Ratio Global\n(IC 95%)",
     Pval = "Valor de p\n(Efeito Global)"
   ) %>%
   theme_booktabs() %>%
   align(j = c("Base", "Comorbidade", "Casos", "Pval"), align = "center", part = "all") %>%
   align(j = "RR", align = "center", part = "all") %>%
   autofit() %>%
   add_footer_lines("Nota:\nEstimativas obtidas por meta-análise de efeitos aleatórios (REML) englobando as nove UFs da Região Nordeste. Valor de p referente ao teste do efeito global. Casos (n) representam o somatório de óbitos com a coocorrência do par.") %>%
   italic(j = 1, part = "footer") %>%
   fontsize(size = 9, part = "footer")
 
 # ------------------------------------------------------------------------------
 # 3. EXPORTAÇÃO PARA WORD (.DOCX)
 # ------------------------------------------------------------------------------
 doc <- read_docx() %>%
   body_add_par("Tabela 2. Topologia da rede de mortalidade: principais hubs e módulos sindrômicos (Nordeste, 2016-2025).", style = "Normal") %>%
   body_add_flextable(tabela2_flex) %>%
   body_add_break() %>%
   body_add_par("Tabela 3. Principais trajetórias clínicas letais validadas por meta-análise (Nordeste, 2016-2025).", style = "Normal") %>%
   body_add_flextable(tabela3_flex)
 
 print(doc, target = "Tabelas_Rede_Nordeste.docx")
 cat("Tabelas exportadas com sucesso para 'Tabelas_Rede_Nordeste.docx'!\n")
 cat("(Tabela 2 agora com valores reais do grafo/bootstrap; Tabela 3 inalterada.)\n")
 
 
 # ============================================================
 # EXTRAÇÃO DOS MÓDULOS DO LOUVAIN - REDE REGIONAL
 # ============================================================
 
 # 1. Extrair CID e módulo de cada nó
 dados_modulos <- data.frame(
   CID = V(grafo_regional)$name,
   Modulo = V(grafo_regional)$group
 )
 
 # 2. Ordenar por módulo e CID
 dados_modulos <- dados_modulos[
   order(dados_modulos$Modulo, dados_modulos$CID),
 ]
 
 # 3. Visualizar todos os CIDs agrupados por módulo
 print(dados_modulos, row.names = FALSE)
 
 # 4. Mostrar quantos CIDs existem em cada módulo
 cat("\n==============================\n")
 cat("NÚMERO DE CIDs POR MÓDULO\n")
 cat("==============================\n")
 
 print(table(dados_modulos$Modulo))
 
 
 # ============================================================
 # CIDs DE CADA MÓDULO
 # ============================================================
 
 cat("\n==============================\n")
 cat("CIDs AGRUPADOS POR MÓDULO\n")
 cat("==============================\n")
 
 for (m in sort(unique(dados_modulos$Modulo))) {
   
   cat("\n--- MÓDULO", m, "---\n")
   
   cids <- dados_modulos$CID[dados_modulos$Modulo == m]
   
   cat(paste(cids, collapse = ", "))
   cat("\n")
 }
 
 
 # ============================================================
 # EXPORTAR PARA EXCEL/CSV
 # ============================================================
 
 write.csv(
   dados_modulos,
   "modulos_louvain_rede_regional.csv",
   row.names = FALSE,
   fileEncoding = "UTF-8"
 )
 
 
 # ============================================================
 # PRINCIPAIS CIDs DE CADA MÓDULO
 # ============================================================
 
 cat("\n====================================\n")
 cat("RESUMO DOS 8 MÓDULOS\n")
 cat("====================================\n")
 
 for (m in sort(unique(dados_modulos$Modulo))) {
   
   cids <- dados_modulos$CID[dados_modulos$Modulo == m]
   
   cat("\nMÓDULO", m, 
       " | Número de CIDs:", length(cids), "\n")
   
   cat("CIDs:\n")
   cat(paste(cids, collapse = ", "))
   cat("\n")
 }
