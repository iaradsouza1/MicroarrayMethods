summarise_probes <- function(df, ref_col, fun = NULL) {

  if(!is.data.frame(df)) {
    stop("Check if 'df' is a data.frame")
  }

  if(!is.null(fun)) {
    id1 <- dplyr::sym(ref_col)
    fun <- dplyr::sym(fun)
    probes <- rownames(df)
    res <- df %>%
      dplyr::mutate(probes = probes) %>%
      tidyr::gather(key = contrast, value = logfc, grep("logfc", colnames(.), ignore.case = T)) %>%
      tidyr::gather(key = degenes, value = de, grep("degene", colnames(.), ignore.case = T)) %>%
      dplyr::group_by((!!id1), contrast) %>%
      dplyr::mutate(logfc = do.call((!!fun), args = list(logfc))) %>%
      dplyr::ungroup() %>%
      dplyr::select((!!id1), grep("probes", colnames(.)),
                    grep("logfc", colnames(.)),
                    grep("degenes", colnames(.)), grep("de", colnames(.))) %>%
      dplyr::filter(de != 0)
    return(res)
  } else {
    id1 <- dplyr::sym(ref_col)
    probes <- rownames(df)
    res <- df %>%
      dplyr::mutate(probes = probes) %>%
      tidyr::gather(key = contrast, value = logfc, grep("logfc", colnames(.), ignore.case = T)) %>%
      tidyr::gather(key = degenes, value = de, grep("degene", colnames(.), ignore.case = T)) %>%
      dplyr::group_by((!!id1), contrast) %>%
      dplyr::slice(which.max(abs(logfc))) %>%
      dplyr::ungroup() %>%
      dplyr::select((!!id1), grep("probes", colnames(.)),
                    grep("logfc", colnames(.)),
                    grep("degenes", colnames(.)), grep("de", colnames(.))) %>%
      dplyr::filter(de != 0)

    return(res)

  }

}


a <- summarise_probes(delimma, "hgnc_symbol", fun =  "mean")
b <- summarise_probes(delimma, "hgnc_symbol", fun = "median")
c <- summarise_probes(delimma, "hgnc_symbol", fun = "var")
d <- summarise_probes(delimma, "hgnc_symbol")


probes <- rownames(delimma)
res <- delimma %>%
  dplyr::mutate(probes = probes) %>%
  tidyr::gather(key = contrast, value = logfc, grep("logfc", colnames(.), ignore.case = T)) %>%
  tidyr::gather(key = degenes, value = de, grep("degene", colnames(.), ignore.case = T)) %>%
  dplyr::group_by(hgnc_symbol, contrast) %>%
  dplyr::slice(which.max(abs(logfc))) %>%
  dplyr::ungroup() %>%
  dplyr::select("hgnc_symbol", grep("probes", colnames(.)),
                grep("logfc", colnames(.)),
                grep("degenes", colnames(.)), grep("de", colnames(.))) %>%
  dplyr::filter(de != 0)






library(dplyr)
library(tidyr)
idx <- grep("log|fc", colnames(delimma), ignore.case = T)
delimma2 <- delimma[,idx]
func <- mean

probes <- rownames(delimma)
a <- delimma %>%
  mutate(probes = probes) %>%
  tidyr::gather(key = contrast, value = logfc, grep("logfc", colnames(delimma), ignore.case = T)) %>%
  group_by(hgnc_symbol, contrast) %>%
  mutate(summary = do.call(func, args = list(logfc))) %>%
  ungroup() %>%
  select(-logfc) %>%
  #mutate(id = 1:n()) %>%
  #spread(contrast, summary) %>%
  select(1, 2, contains("summary"), contains("degenes"), contains("probes"))

a <- delimma %>%
  dplyr::mutate(probes = probes) %>%
  tidyr::gather(key = contrast, value = logfc, grep("logfc", colnames(.), ignore.case = T)) %>%
  tidyr::gather(key = degenes, value = de, grep("degene", colnames(.), ignore.case = T)) %>%
  dplyr::group_by(hgnc_symbol, contrast) %>%
  dplyr::mutate(logfc = do.call(func, args = list(logfc))) %>%
  dplyr::ungroup() %>%
  #mutate(id = 1:n()) %>%
  #spread(contrast, summary) %>%
  dplyr::select(1,2, contains("probes"), contains("logFC"), contains("degenes"), contains("de")) %>%
  filter(de != 0)


a <- delimma %>%
  gather(key = contrast, value = logfc, grep("log|fc", colnames(delimma), ignore.case = T)) %>%
  group_by(hgnc_symbol, contrast) %>%
  slice(which.max(abs(logfc))) %>%
  ungroup() %>%
  spread(contrast, logfc) %>%
  select(1,2,contains("logFC"), contains("degenes"))

