Topologia da Mortalidade Hospitalar no Nordeste (2016–2025)

Este repositório contém o pipeline analítico, topológico e de ciência aberta para o estudo da mortalidade hospitalar na Região Nordeste do Brasil (2016–2025).

O código extrai dados reais do SIH/SUS, purifica as coortes clínicas e constrói redes de associações de comorbidades letais utilizando meta-análise, teoria dos grafos e simulações de Monte Carlo.

📊 Sobre o Projeto

O script em R fornecido neste repositório executa um pipeline completo (Versão 4.0/4.1) dividido em módulos sequenciais, incluindo:

Extração e Purificação: Download automatizado de microdados do SIH-RD (DATASUS) para as 9 UFs do Nordeste.

Meta-Análise entre Estados: Avaliação de associações clínicas (Diagnóstico Base x Comorbidade) com modelos de efeitos aleatórios (REML/Knapp-Hartung).

Topologia de Redes e Modelos Nulos: Comparação da rede observada com modelos Erdős–Rényi, Barabási–Albert e Watts–Strogatz.

Simulação de Monte Carlo: Modelo estocástico de "cascata sindrômica" (percolação de aresta) para identificar diagnósticos-fonte de alto impacto.

Análise de Robustez e Subgrupos: Testes de resiliência da rede, impacto das eras epidemiológicas da COVID-19 (Pré, Auge e Pós-pandemia) e análises de sensibilidade.

⚙️ Pré-requisitos e Instalação

Para reproduzir este estudo, você precisará do R instalado. O script automaticamente checa e instala a maioria das dependências, mas certifique-se de ter os seguintes pacotes principais disponíveis:

install.packages(c("tidyverse", "igraph", "metafor", "lme4", "poweRlaw", "visNetwork", "ggraph", "patchwork", "DiagrammeR", "geobr", "sf"))


Nota: O pacote microdatasus será instalado automaticamente via GitHub caso não esteja presente.

🚀 Como Executar

Clone o repositório:

git clone https://github.com/SEU-USUARIO/NOME-DO-REPOSITORIO.git


Abra o script principal (pipeline_mortalidade_nordeste.R ou o nome que você definiu para o seu .R).

Execute o código sequencialmente.

⚠️ Sistema de Checkpoints (Importante)

O script processa milhões de registros do DATASUS. Para evitar perda de progresso, o código possui um sistema de checkpointing.

Ao final de cada módulo pesado (download, purificação, cálculo 2x2, meta-análise), os resultados são salvos em uma pasta local chamada checkpoints_nordeste/.

Se o script travar ou for interrompido, basta rodá-lo novamente. Ele pulará as etapas já concluídas lendo os arquivos do disco.

Para forçar o reprocessamento de uma etapa, apague o arquivo .rds correspondente na pasta de checkpoints.

📁 Estrutura de Saída (Outputs)

Após a execução completa, o pipeline irá gerar na raiz do projeto:

Figuras (PNG/TIFF): Mapas de calor (Cartografia do NE), Forest Plots (estilo Cochrane), Histogramas de Monte Carlo, Fluxograma STROBE e Grafos Topológicos.

Tabelas (.docx / .csv): Tabelas demográficas prontas para publicação (formato APA) e métricas descritivas.

Workbook OSF (.xlsx): Uma planilha completa (OSF_Pacote_Completo_Topologia_Mortalidade_NORDESTE_v4.1.xlsx) contendo abas com todas as métricas detalhadas (I², Risco Relativo, Centralidades, Sensibilidade, etc.) pronta para submissão a repositórios de ciência aberta.
