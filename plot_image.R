plot_image<- function(image_vector,ncol) { 
  plot(as.cimg(t(matrix(image_vector, ncol=ncol))), axes=FALSE, asp=1,xlab=xlab)
}
