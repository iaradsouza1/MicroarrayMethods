
markdown <- c(
  '---',
    'title: "Quality control"',
  'output: ',
    '\t html_document: ',
    '\t\t theme: spacelab',
  '\t highlight: tango',
  '---',

  '# Raw data',

  '```{r, echo=FALSE, message=FALSE, warning=FALSE}',
  'raw <- ReadAffy(celfile.path = "~/Área de Trabalho/Cassia_3/Analise_2/GSE5617_RAW/")',
  'eset_raw <- ExpressionSet(assayData = raw@assayData$exprs)',
  '```',

  '```{r, echo=FALSE, message=FALSE, warning=FALSE}',
  'plotDensities(eset_raw, legend = FALSE, main = "Raw data samples density plot")',
  '```',

  '```{r, echo=FALSE, message=FALSE, warning=FALSE}',
  'boxplot(raw, main = "Samples boxplot", las = 2)',
  '```',

  '```{r, echo=FALSE, message=FALSE, warning=FALSE}',
  'par(mfrow = c(1,1))',
  'deg <- AffyRNAdeg(raw)',
  'plotAffyRNAdeg(deg)',
  '```',

  '# Normalized data',

  '```{r, echo=FALSE, message=FALSE, warning=FALSE}',
  'eset_normalized <- rma(raw)',
  'affy::MAplot(eset_raw, plot.method = "smoothScatter")',
  '```'



)

markdown::markdownToHTML(text = knitr::knit(text = markdown), output = 'teste2.html')

