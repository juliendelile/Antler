.onLoad <- function(libname, pkgname) {
  op <- options()
  set.seed(0)
  op.antler <- list(
    antler.feature.colors      = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]),
    antler.cell_cluster.colors = c(
    	RColorBrewer::brewer.pal(9, "Set1"),
    	RColorBrewer::brewer.pal(12, "Set3"),
    	'#EBFF00','#8F00FF','#14FF00','#00FFFF','#24FF00','#FF00E6','#FF000F','#FF3D00',
    	'#3300FF','#00B2FF','#1400FF','#0075FF','#FF7A00','#9E00FF','#FF8A00','#CCFF00',
    	'#0066FF','#FF0F00','#FF9900','#05FF00','#FF004C','#FFF500','#FF005C','#0047FF',
    	'#FFD600','#00A3FF','#2400FF','#FF00D6','#00FF57','#00FF47','#FF001F','#FF00B8',
    	'#ADFF00','#00FF75','#70FF00','#FF0099','#00FFA3','#00FF29','#FFE500','#0500FF',
    	'#FF6B00','#00E0FF','#00FF1A','#0085FF','#00FFF0','#FF007A','#BD00FF','#0057FF',
    	'#8000FF','#00FFE0','#00FF85','#FF00C7','#FFA800','#0094FF','#FFC700','#DBFF00',
    	'#EB00FF','#00FFC2','#BDFF00','#33FF00','#0038FF','#FF0000','#FA00FF','#FF003D',
    	'#00C2FF','#4200FF','#52FF00','#00FFD1','#9EFF00','#00FF94','#8FFF00','#FF5C00',
    	'#FFB800','#FF002E','#000AFF','#FF00A8','#CC00FF','#00FF38','#00FF66','#00FFB3',
    	'#00FF0A','#FF00F5','#00D1FF','#00F0FF','#5200FF','#0019FF','#80FF00','#7000FF',
    	'#0029FF','#FF1F00','#FF008A','#FF4D00','#AD00FF','#FAFF00','#DB00FF','#42FF00',
    	'#6100FF','#FF006B','#FF2E00','#61FF00')
  )
  toset <- !(names(op.antler) %in% names(op))
  if(any(toset)) options(op.antler[toset])

  invisible()
}