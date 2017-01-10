
# This file contains the code to load the libraries and functions required by MAMO and associated functions

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# load libraries

library(truncnorm)
library(compiler)
library(pastecs)
library(lme4)
library(plotrix)
library(rsm)
library(ade4)
library(ggplot2)
library(msm)

# below are grouped all functions of the package factoextra
# It was used here for plotting correlation circles of PCA; the package proved difficult to obtain via 'classical' R routes

#' Extract and visualize the eigenvalues/variances of dimensions
#' 
#' @description
#' Extracts and plots the eigenvalues/variances of the dimensions 
#' from the results of Principal Component Analysis (PCA), 
#' Correspondence Analysis (CA) and 
#' Multiple Correspondence Analysis (MCA) functions.\cr\cr
#' \itemize{
#' \item{get_eig(): Extract the eigenvalues/variances of the principal dimensions}
#' \item{fviz_eig(): Plot the eigenvalues/variances against the number of dimensions}
#' \item{get_eigenvalue(): an alias of get_eig()}
#' \item{fviz_screeplot(): an alias of fviz_eig()}
#' }
#' 
#' @param X an object of class PCA, CA and MCA [FactoMineR]; prcomp and princomp [stats]; 
#'  dudi, pca, coa and acm [ade4]; ca and mjca [ca package].
#' @param choice a text specifying the data to be plotted. 
#' Allowed values are "variance" or "eigenvalue".
#' @param geom a text specifying the geometry to be used for the graph.
#'  Allowed values are "bar" for barplot, "line" for lineplot or c("bar", "line") to use both types.
#' @param barfill fill color for bar plot.
#' @param barcolor outline color for bar plot.
#' @param linecolor color for line plot (when geom contains "line").
#' @param ncp a numeric value specifying the number of dimensions to be shown.
#' @param addlabels logical value. If TRUE, labels are added at the top of bars or points
#'  showing the information retained by each dimension.
#'  @param ... optional arguments to be passed to the functions geom_bar(), 
#'  geom_line(), geom_text() or fviz_eig().
#'  
#' @return 
#' \itemize{
#' \item{get_eig() (or get_eigenvalue()): returns a data.frame containing 3 columns: 
#' the eigenvalues, the percentage of variance and  the cumulative percentage of variance 
#' retained by each dimension.}
#' \item{fviz_eig() (or fviz_screeplot()): returns a ggplot2}
#' }
#'  
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal Component Analysis
#' # ++++++++++++++++++++++++++
#' data(iris)
#' res.pca <- prcomp(iris[, -5],  scale = TRUE)
#' 
#' # Extract eigenvalues/variances
#' get_eig(res.pca)
#' 
#' # Default plot
#' fviz_eig(res.pca)
#' 
#' # Add labels
#' fviz_eig(res.pca, addlabels=TRUE)
#' 
#' # Change the y axis limits
#' fviz_eig(res.pca, addlabels=TRUE, hjust = -0.3) +
#'    ylim(0, 80)
#' # Scree plot - Eigenvalues
#' fviz_eig(res.pca, choice = "eigenvalue", addlabels=TRUE)
#' 
#' # Use only bar plot
#' fviz_eig(res.pca, geom="bar", width=0.8)
#' 
#' # Use only line plot
#' fviz_eig(res.pca, geom="line")
#' 
#' # Change theme
#' fviz_screeplot(res.pca) + theme_minimal()
#' # theme_classic()
#' fviz_eig(res.pca) + theme_classic()
#' 
#' # Customized plot
#' fviz_eig(res.pca, addlabels=TRUE, hjust = -0.3,
#'            linecolor ="red") + theme_minimal()
#' # Change colors, y axis limits and theme           
#' p <- fviz_eig(res.pca, addlabels=TRUE, hjust = -0.3,
#'                barfill="white", barcolor ="darkblue",
#'                linecolor ="red") + ylim(0, 85) + 
#'                theme_minimal()
#' print(p)
#' # Change titles
#' p + labs(title = "Variances - PCA",
#'         x = "Principal Components", y = "% of variances")
#'         
#' # Correspondence Analysis
#' # +++++++++++++++++++++++++++++++++
#' library(FactoMineR)
#' data(housetasks)
#' res.ca <- CA(housetasks, graph = FALSE)
#' get_eig(res.ca)
#' fviz_eig(res.ca)
#' 
#' # Multiple Correspondence Analysis
#' # +++++++++++++++++++++++++++++++++
#' library(FactoMineR)
#' data(poison)
#' res.mca <- MCA(poison, quanti.sup = 1:2, 
#'               quali.sup = 3:4, graph=FALSE)
#' get_eig(res.mca)
#' fviz_eig(res.mca)
#'  }
#'
#' @name eigenvalue
NULL
#' @rdname eigenvalue
#' @export
get_eig<-function(X){
  
  # FactoMineR package
  if(inherits(X, c('PCA', 'CA', 'MCA'))) eig <- X$eig
  else{
    # stats package
    if(inherits(X, 'prcomp') | inherits(X, 'princomp')) eig <- (X$sdev)^2
    # ade4 package
    else if(inherits(X, c('pca', 'coa', 'acm')) & inherits(X, 'dudi')) eig <- X$eig
    # ca package
    else if(inherits(X, 'ca'))  eig <- X$sv^2
    else if(inherits(X, 'mjca')) eig <- X$inertia.e
    # MASS
    else if(inherits(X, 'correspondence'))  eig <- X$cor^2
    else stop("An object of class : ", class(X), 
              " can't be handled by the function get_eigenvalue()")
    
    variance <- eig*100/sum(eig)
    cumvar <- cumsum(variance)
    eig <- data.frame(eigenvalue = eig, variance = variance, 
                      cumvariance = cumvar)
  }
  
  colnames(eig) <- c("eigenvalue", "variance.percent", 
                     "cumulative.variance.percent")
  rownames(eig) <- paste0("Dim.", 1:nrow(eig))
  
  eig 
}

#' @rdname eigenvalue
#' @export
get_eigenvalue <- function(X){
  get_eig(X)
}

#' @rdname eigenvalue
#' @export
fviz_eig<-function(X, choice=c("variance", "eigenvalue"), geom=c("bar", "line"),
                         barfill="steelblue", barcolor="steelblue", linecolor = "black",
                         ncp=10, addlabels=FALSE, ...)
{
  
  eig <- get_eigenvalue(X)
  eig <-eig[1:min(ncp, nrow(eig)), , drop=FALSE]
  
  title <- "Scree plot"
  xlab <- "Dimensions"
  ylab <- "Percentage of explained variances"
  
  choice <- choice[1]
  if(choice=="eigenvalue") {
    eig <- eig[,1]
    text_labels <- round(eig,1)
    ylab <- "Eigenvalue"
  }
  else if(choice=="variance") {
    eig <- eig[,2]
    text_labels <- paste0(round(eig,1), "%")
  }
  else stop("Allowed values for the argument choice are : 'variance' or 'eigenvalue'")
  
  if(length(intersect(geom, c("bar", "line"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  
  df.eig <- data.frame(dim = factor(1:length(eig)), eig=eig )
  p <- ggplot(df.eig, aes(dim, eig, group=1 ))
  if("bar" %in% geom) p <- p + geom_bar(stat="identity", fill=barfill, color = barcolor,...)
  if("line" %in% geom) p <- p + geom_line(color = linecolor, ...)+
    geom_point(shape=19, color=linecolor)
  if(addlabels) p <- p + geom_text(label = text_labels,
                                   vjust=-0.4, ...)
  p <- p + labs(title = title, x = xlab, y = ylab)
  
  p 
}


#' @rdname eigenvalue
#' @export 
fviz_screeplot<- function(...){
  fviz_eig(...)
} 

#' @include utilities.R get_pca.R eigenvalue.R
NULL
#' Subset and summarize the output of factor analyses
#' 
#' @description
#'  Subset and summarize the results of Principal Component Analysis (PCA), 
#' Correspondence Analysis (CA) and 
#' Multiple Correspondence Analysis (MCA) functions from several packages.
#' @param X an object of class PCA, CA and MCA [FactoMineR]; prcomp and princomp [stats]; 
#'  dudi, pca, coa and acm [ade4]; ca [ca package].
#' @param element the element to subset from the output. Possible values are
#'  "row" or "col" for CA; "var" or "ind" for PCA and MCA
#' @param result the result to be extracted for the element. Possible values are
#'  the combination of c("cos2", "contrib", "coord")
#' @param axes a numeric vector specifying the axes of interest. Default values are 1:2
#'  for axes 1 and 2.
#' @param select a selection of variables. Allowed values are NULL or a list containing the arguments
#'  name, cos2 or contrib. Default is list(name = NULL, cos2 = NULL, contrib = NULL):
#'  \itemize{
#'  \item name: is a character vector containing variable names to be selected
#'  \item cos2: if cos2 is in [0, 1], ex: 0.6, then variables with a cos2 > 0.6 are selected.
#'   if cos2 > 1, ex: 5, then the top 5 variables with the highest cos2 are selected
#' \item contrib: if contrib > 1, ex: 5,  then the top 5 variables with the highest cos2 are selected. 
#'  }
#' @return A data frame containing the (total) coord, cos2 and the contribution for the axes.
#' @details If length(axes) > 1, then the columns contrib and cos2 correspond to the total contributions and total cos2
#'  of the axes. In this case, the column coord is calculated as x^2 + y^2 + ...+; x, y, ... are the coordinates of
#'  the points on the specified axes.
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal component analysis
#' # +++++++++++++++++++++++++++++
#' data(decathlon2)
#' decathlon2.active <- decathlon2[1:23, 1:10]
#' res.pca <- prcomp(decathlon2.active,  scale = TRUE)
#' 
#' # Summarize variables on axes 1:2
#' facto_summarize(res.pca, "var", axes = 1:2)[,-1]
#' # Select the top 5 contributing variables
#' facto_summarize(res.pca, "var", axes = 1:2,
#'            select = list(contrib = 5))[,-1]
#' # Select variables with cos2 >= 0.6
#' facto_summarize(res.pca, "var", axes = 1:2,
#'            select = list(cos2 = 0.6))[,-1]
#' # Select by names
#' facto_summarize(res.pca, "var", axes = 1:2,
#'      select = list(name = c("X100m", "Discus", "Javeline")))[,-1]
#'            
#' # Summarize individuals on axes 1:2
#' facto_summarize(res.pca, "ind", axes = 1:2)[,-1]
#' 
#' # Correspondence Analysis
#' # ++++++++++++++++++++++++++
#' # Install and load FactoMineR to compute CA
#' # install.packages("FactoMineR")
#' library("FactoMineR")
#' data("housetasks")
#' res.ca <- CA(housetasks, graph = FALSE)
#' # Summarize row variables on axes 1:2
#' facto_summarize(res.ca, "row", axes = 1:2)[,-1]
#' # Summarize column variables on axes 1:2
#' facto_summarize(res.ca, "col", axes = 1:2)[,-1]
#' 
#' # Multiple Correspondence Analysis
#' # +++++++++++++++++++++++++++++++++
#' library(FactoMineR)
#' data(poison)
#' res.mca <- MCA(poison, quanti.sup = 1:2, 
#'               quali.sup = 3:4, graph=FALSE)
#' # Summarize variables on axes 1:2
#' res <- facto_summarize(res.mca, "var", axes = 1:2)
#' head(res)
#' # Summarize individuals on axes 1:2
#' res <- facto_summarize(res.mca, "ind", axes = 1:2)
#' head(res)
#' 
#'  }
#' @export 
facto_summarize <- function(X, element,
                            result = c("coord", "cos2", "contrib"),
                            axes=1:2, select = NULL)
                            
  { 
  # check element
  if(!element %in% c("row", "col", "var", "ind"))
    stop('Te argument element should be one of "row", "col", "var", "ind"')
  
  # check and get the classe of X
  facto_class <- .get_facto_class(X)
  
  # Extract the element
  element <- element[1]
  if(facto_class=="CA"){
  if(element %in% c("ind", "row")) elmt<- get_ca_row(X)
  else if(element  %in% c("var", "col") ) elmt <- get_ca_col(X)
  }
  else if(facto_class=="PCA"){
    if(element %in% c("var", "col")) elmt<- get_pca_var(X)
    else if(element %in% c("ind", "row")) elmt <- get_pca_ind(X)
  }
  else if(facto_class=="MCA"){
    if(element %in% c("var", "col")) elmt<- get_mca_var(X)
    else if(element %in% c("ind", "row")) elmt <- get_mca_ind(X)
  }
  
  
  # check axes
  if(max(axes) > ncol(elmt$coord))
    stop("The value of the argument axes is incorrect. ",
         "The number of axes in the data is: ", ncol(elmt$coord), 
         ". Please try again with axes between 1 - ", ncol(elmt$coord))
  
  # summarize the result
  res = NULL
  
  # 1.Extract the coordinates x, y and coord
  if("coord" %in% result){
    dd <- data.frame(elmt$coord[, axes, drop=FALSE])
    coord <- apply(dd^2, 1, sum) # x^2 + y2 + ...
    res = cbind(dd, coord = coord)
  }
  
  # 2. Extract the cos2
  if("cos2" %in% result){
    cos2 <- elmt$cos2[, axes]
    if(length(axes) > 1) cos2 <- apply(cos2, 1, sum, na.rm=TRUE)
    res <- cbind(res, cos2 = cos2)
  }
  
  # 3. Extract the contribution
  if("contrib" %in% result){
    contrib <- elmt$contrib[, axes]
    if(length(axes) > 1) {
      eig <- get_eigenvalue(X)[axes,1]
      # Adjust variable contributions by the Dimension eigenvalues
      contrib <- t(apply(contrib, 1, 
                         function(var.contrib, pc.eig){var.contrib*pc.eig},
                         eig))
      contrib <-apply(contrib, 1, sum)
    }
    res <- cbind(res, contrib = contrib)
  }
  
  name <- rownames(elmt$coord)
  if(is.null(name)) name <- as.character(1:nrow(elmt$coord))
  res <- cbind.data.frame(name = name, res)
  if(!is.null(select)) res <- .select(res, select)
  res 
}

#' Add supplementary data to a plot
#' 
#' @description
#' Add supplementary data to a plot
#'  
#' @param ggp a ggplot2 plot.
#' @param df a data frame containing the x and y coordinates
#' @param axes a numeric vector of length 2 specifying the components to be plotted.
#' @param geom a character specifying the geometry to be used for the graph
#'  Allowed values are "point" or "arrow" or "text"
#' @param color the color to be used
#' @param addlabel a logical value. If TRUE, labels are added
#' @param labelsize the size of labels. Default value is 4
#' @param pointsize the size of points
#' @param shape point shape when geom ="point"
#' @param linetype the linetype to be used when geom ="arrow"
#' @param jitter a parameter used to jitter the points in order to reduce overplotting. 
#' It's a list containing the objects what, width and height (i.e jitter = list(what, width, height)). 
#' \itemize{
#' \item what: the element to be jittered. Possible values are "point" or "p"; "label" or "l"; "both" or "b".
#' \item width: degree of jitter in x direction
#' \item height: degree of jitter in y direction
#' }
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal component analysis
#' data(decathlon2)
#' decathlon2.active <- decathlon2[1:23, 1:10]
#' res.pca <- prcomp(decathlon2.active,  scale = TRUE)
#' 
#' # Visualize variables
#' p <- fviz_pca_var(res.pca)
#' print(p)
#' 
#' # Add supplementary variables
#' coord <- data.frame(PC1 = c(-0.7, 0.9), PC2 = c(0.25, -0.07))
#' rownames(coord) <- c("Rank", "Points")
#' print(coord)
#' fviz_add(p, coord, color ="blue", geom="arrow") 
#'  }
#'  
#' @export 
fviz_add <- function(ggp, df, axes = c(1,2), geom=c("point", "arrow"), color ="blue", 
                     addlabel = TRUE, labelsize = 4, pointsize = 2, shape=19, linetype ="dashed",
                     jitter = list(what = "label", width = NULL, height = NULL))
{
  if(!inherits(df, c("data.frame", "matrix")))
     stop("df should be a data frame or a matrix")
     
  if(ncol(df) < 2)
    stop("df should have at least two columns (x and y coordinates)")
  
  if(length(intersect(geom, c("point", "arrow", "text"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  
  df <- data.frame(name = rownames(df), x = df[,axes[1]], y = df[,axes[2]]) 
  label_coord <- df
  
  # jittering
  if(jitter$what %in% c("both", "b")){
    label_coord <- df <- .jitter(data, jitter)
  }
  else if(jitter$what %in% c("point", "p")){
    df<- .jitter(df, jitter)
  }
  else if(jitter$what %in% c("label", "l")){
    label_coord <- .jitter(label_coord, jitter)
  }
  
  
  if("point" %in% geom) {
    p <-  ggp + geom_point(data = df, aes_string("x", "y"), 
                           color = color, shape = shape, size = pointsize)
    if(addlabel) 
      p <- p + geom_text(data = label_coord, aes_string("x", "y"), color = color,
                         label = df$name, size = labelsize, vjust=-0.7)
  }
  else if("arrow" %in% geom){
    p <- ggp + geom_segment(data = df,
                      aes_string(x = 0, y = 0, xend = 'x', yend = 'y'),
                      arrow = grid::arrow(length = grid::unit(0.2, 'cm')), 
                      color=color, linetype=linetype)
    if(addlabel)
      p <- p + geom_text(data = label_coord, aes_string("x", "y"),
                         label = df$name, color=color, 
                         size = labelsize, hjust=0.8, vjust=0) 
  }
  else if("text" %in% geom)
    p <- ggp + geom_text(data = label_coord, aes_string("x", "y"), color = color,
                       label = df$name, size = labelsize, vjust=-0.7)
  
  return(p)
}

#' @include utilities.R
NULL
#' Visualize Correspondence Analysis
#' 
#' @description
#' Graph of column/row variables from the output of Correspondence Analysis (CA).\cr\cr
#' \itemize{
#' \item{fviz_ca_row(): Graph of row variables}
#' \item{fviz_ca_col(): Graph of column variables}
#' \item{fviz_ca_biplot(): Biplot of row and column variables}
#' \item{fviz_ca(): An alias of fviz_ca_biplot()}
#' } 
#' @param X an object of class CA [FactoMineR], ca [ca], coa [ade4];
#'  correspondence [MASS].
#' @param axes a numeric vector of length 2 specifying the dimensions to be plotted.
#' @param shape.row,shape.col the point shapes to be used for row/column variables. 
#' Default values are 19 for rows and 17 for columns.
#' @param geom a character specifying the geometry to be used for the graph. 
#' Allowed values are the combination of c("point", "arrow", "text"). 
#' Use "point" (to show only points); 
#' "text" to show only labels; c("point", "text") or c("arrow", "text") to show both types.
#' @param label a character vector specifying the elements to be labelled. 
#' Default value is "all". Allowed values are "none" or the combination of c("row", "row.sup", "col", "col.sup"). 
#' Use "col" to label only active column variables; "col.sup" to label only supplementary columns; etc
#' @param invisible a character value specifying the elements to be hidden on the plot. 
#' Default value is "none". Allowed values are 
#' the combination of c("row", "row.sup","col", "col.sup").
#' @param labelsize font size for the labels
#' @param pointsize the size of points
#' @param col.col,col.row color for column/row points. 
#' The default values are "red" and "blue", respectively. 
#' Allowed values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the colors for row/column variables are automatically controlled by their qualities ("cos2"),
#'  contributions ("contrib"), coordinates (x^2 + y^2, "coord"), x values("x") or y values("y")
#' @param alpha.col,alpha.row controls the transparency of colors.
#' The value can variate from 0 (total transparency) to 1 (no transparency).
#' Default value is 1. Allowed values include also : "cos2", "contrib", "coord", "x" or "y" 
#' as for the arguments col.col and col.row.
#' @param col.col.sup,col.row.sup colors for the supplementary column and row points, respectively.
#' @param select.col,select.row a selection of columns/rows to be drawn. 
#' Allowed values are NULL or a list containing the arguments name, cos2 or contrib: 
#' \itemize{
#' \item name is a character vector containing column/row names to be drawn
#' \item cos2 if cos2 is in [0, 1], ex: 0.6, then columns/rows with a cos2 > 0.6 are drawn. 
#' if cos2 > 1, ex: 5, then the top 5 columns/rows with the highest cos2 are drawn.
#' \item contrib if contrib > 1, ex: 5,  then the top 5 columns/rows with the highest cos2 are drawn
#' }
#' @param map character string specifying the map type. Allowed options include: 
#' "symmetric", "rowprincipal", "colprincipal", "symbiplot", "rowgab", 
#' "colgab", "rowgreen" and "colgreen". See details
#' @param arrows Vector of two logicals specifying if the plot should contain
#'  points (FALSE, default) or arrows (TRUE).
#'  First value sets the rows and the second value sets the columns.
#' @param jitter a parameter used to jitter the points in order to reduce overplotting. 
#' It's a list containing the objects what, width and height (i.e jitter = list(what, width, height)). 
#' \itemize{
#' \item what: the element to be jittered. Possible values are "point" or "p"; "label" or "l"; "both" or "b".
#' \item width: degree of jitter in x direction
#' \item height: degree of jitter in y direction
#' }
#' @param ... optional arguments.
#' @details The default plot of CA is a "symmetric" plot in which both rows and 
#' columns are in principal coordinates. In this situation, it's not possible 
#' to interpret the distance between row points and column points. To overcome this 
#' problem, the simplest way is to make an asymmetric plot. This means that, 
#' the column profiles must be presented in row space or vice-versa. 
#' The allowed options for the argument map are: 
#' \itemize{
#' \item "rowprincipal" or "colprincipal": asymmetric plots with either rows in principal 
#' coordinates and columns in standard coordinates, or vice versa. 
#' These plots preserve row metric or column metric respectively.
#' \item "symbiplot": Both rows and columns are scaled to have variances 
#' equal to the singular values (square roots of eigenvalues), 
#' which gives a symmetric biplot but does not preserve row or column metrics. 
#' \item "rowgab" or "colgab": Asymmetric maps, proposed by Gabriel & Odoroff (1990), 
#' with rows (respectively, columns) in 
#' principal coordinates and columns (respectively, rows) in standard coordinates 
#' multiplied by the mass of the corresponding point.
#' \item "rowgreen" or "colgreen": The so-called contribution biplots 
#' showing visually the most contributing points (Greenacre 2006b). 
#' These are similar to "rowgab" and "colgab" except that the points 
#' in standard coordinates are multiplied by the square root of the corresponding masses, 
#' giving reconstructions of the standardized residuals.
#' }
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Correspondence Analysis
#' # ++++++++++++++++++++++++++++++
#' # Install and load FactoMineR to compute CA
#' # install.packages("FactoMineR")
#' library("FactoMineR")
#' data(housetasks)
#' head(housetasks)
#' res.ca <- CA(housetasks, graph=FALSE)
#' 
#' # Graph of row variables
#' # +++++++++++++++++++++
#' # Default plot
#' fviz_ca_row(res.ca)
#' # Change title and axis labels
#' fviz_ca_row(res.ca) +
#'  labs(title = "CA", x = "Dim.1", y ="Dim.2" )
#' # Change axis limits by specifying the min and max
#' fviz_ca_row(res.ca) + 
#'    xlim(-1.3, 1.7) + ylim (-1.5, 1)
#' # Use text only
#' fviz_ca_row(res.ca, geom = "text")
#' # Use points only
#' fviz_ca_row(res.ca, geom="point")
#' # Change the size of points
#' fviz_ca_row(res.ca, geom="point", pointsize = 4)
#' # Change point color and theme
#' fviz_ca_row(res.ca, col.row = "violet")+
#'    theme_minimal()
#'    
#' # Control automatically the color of row points
#' # using the cos2 or the contributions
#' # cos2 = the quality of the rows on the factor map
#' fviz_ca_row(res.ca, col.row="cos2") 
#' # Gradient color
#' fviz_ca_row(res.ca, col.row="cos2") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=0.5)
#' # Change the theme and use only points
#' fviz_ca_row(res.ca, col.row="cos2", geom = "point") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=0.4)+ theme_minimal()
#'       
#' # Color by the contributions   
#' fviz_ca_row(res.ca, col.row="contrib") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=10)
#'       
#' # Control the transparency of the color by the
#' # contributions
#' fviz_ca_row(res.ca, alpha.row="contrib") +
#'      theme_minimal()        
#'              
#' # Select and visualize rows with cos2 > 0.5
#' fviz_ca_row(res.ca, select.row = list(cos2 = 0.5))
#' # Select the top 7 according to the cos2
#' fviz_ca_row(res.ca, select.row = list(cos2 = 7))
#' # Select the top 7 contributing rows
#' fviz_ca_row(res.ca, select.row = list(contrib = 7))
#' # Select by names
#' fviz_ca_row(res.ca, 
#' select.row = list(name = c("Breakfeast", "Repairs", "Holidays")))
#' 
#'  
#' # Graph of column points
#' # ++++++++++++++++++++++++++++
#' # Default plot
#' fviz_ca_col(res.ca)
#' # Change color and theme
#' fviz_ca_col(res.ca, col.col="steelblue")+
#'  theme_minimal()
#'  
#' # Control colors using their contributions
#' fviz_ca_col(res.ca, col.col = "contrib")+
#'  scale_color_gradient2(low = "white", mid = "blue", 
#'            high = "red", midpoint = 25) +
#'  theme_minimal()          
#' # Control the transparency of variables using their contributions
#' fviz_ca_col(res.ca, alpha.col = "contrib") +
#'    theme_minimal()
#'    
#' # Select and visualize columns with cos2 >= 0.4
#' fviz_ca_col(res.ca, select.col = list(cos2 = 0.4))
#' # Select the top 3 contributing columns
#' fviz_ca_col(res.ca, select.col = list(contrib = 3))
#' # Select by names
#' fviz_ca_col(res.ca, 
#'  select.col= list(name = c("Wife", "Husband", "Jointly")))
#'     
#' # biplot
#' # ++++++++++++++++++++++++++
#' # Symetric Biplot of rows and columns
#' fviz_ca_biplot(res.ca)
#' # Asymetric biplot, use arrows for columns
#' fviz_ca_biplot(res.ca, map ="rowprincipal",
#'  arrow = c(FALSE, TRUE))
#' # Keep only the labels for row points
#' fviz_ca_biplot(res.ca, label ="row")
#' # Keep only labels for column points
#' fviz_ca_biplot(res.ca, label ="col")
#' # Hide row points
#' fviz_ca_biplot(res.ca, invisible ="row")
#' # Hide column points
#' fviz_ca_biplot(res.ca, invisible ="col")
#'# Control automatically the color of rows using the cos2
#' fviz_ca_biplot(res.ca, col.row="cos2") +
#'        theme_minimal()
#' # Select the top 7 contributing rows
#' # And the top 3 columns
#' fviz_ca_biplot(res.ca,  
#'                select.row = list(contrib = 7),
#'                select.col = list(contrib = 3)) 
#' }
#'  
#' @name fviz_ca
#' 
#' @rdname fviz_ca 
#' @export 
fviz_ca_row <-function(X,  axes = c(1,2), shape.row = 19, 
                       geom=c("point", "text"),
                       label = "all", invisible="none", labelsize=4, pointsize = 2,
                       col.row ="blue", col.row.sup="darkblue",  alpha.row = 1,
                       select.row = list(name = NULL, cos2 = NULL, contrib = NULL),
                       map ="symmetric",
                       jitter = list(what = "label", width = NULL, height = NULL),...)
{
  
  if(length(intersect(geom, c("point", "text", "arrow"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  if(length(axes) > 2) stop("axes should be of length 2")
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  row <- facto_summarize(X, element = "row", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(row)[2:3] <-  c("x", "y")
  
  # scale row coords according to the type of map
  row <- .scale_ca(row, res.ca = X,  element = "row", 
                   type = map, axes = axes)
  
  row.all <- row
  if(!is.null(select.row)) row <- .select(row, select.row)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.row %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(row.all[, alpha.row])
  
  p <- ggplot() 
  if(hide$row) p <-ggplot()+geom_blank(data=row, aes_string("x","y"))
  else p <- .ggscatter(data = row, x = 'x', y = 'y', 
                       col=col.row,  alpha = alpha.row, 
                       alpha.limits = alpha.limits, shape = shape.row, 
                       geom = geom, lab = lab$row, labelsize = labelsize,
                       pointsize = pointsize, jitter = jitter)
  
  # Add supplementary rows
  if(inherits(X, c('CA', 'ca')) & !hide$row.sup){
    row_sup <- .get_supp(X, element = "row.sup", axes = axes,
                         select = select.row)
    if(!is.null(row_sup)){
      colnames(row_sup)[2:3] <-  c("x", "y")
      row_sup <- .scale_ca(row_sup, res.ca = X,  element = "row.sup", 
                           type = map, axes = axes)
    }
    
    if(!is.null(row_sup)){
      p <- fviz_add(p, df = row_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.row.sup, shape = shape.row,
                    labelsize = labelsize, addlabel = (lab$row.sup & "text" %in% geom),
                    pointsize = pointsize, jitter = jitter)
      
    }   
  } 
  
  p <- .fviz_finish(p, X, axes)
  p
  
}

#' @rdname fviz_ca
#' @export 
fviz_ca_col <-function(X,  axes = c(1,2), shape.col = 17, 
                       geom=c("point", "text"),
                       label = "all", invisible="none", labelsize=4, pointsize = 2,
                       col.col ="red", col.col.sup="darkred",  alpha.col = 1,
                       select.col = list(name = NULL, cos2 = NULL, contrib = NULL),
                       map ="symmetric",
                       jitter = list(what = "label", width = NULL, height = NULL), ...)
{
  
  if(length(intersect(geom, c("point", "text", "arrow"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed")
  if(length(axes) > 2) stop("axes should be of length 2")
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  col <- facto_summarize(X, element = "col", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(col)[2:3] <-  c("x", "y")
  
  # scale coords according to the type of map
  col <- .scale_ca(col, res.ca = X,  element = "col", 
                   type = map, axes = axes)
  
  col.all <- col
  if(!is.null(select.col)) col <- .select(col, select.col)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.col %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(col.all[, alpha.col])
  
  p <- ggplot() 
  if(hide$col) p <-ggplot()+geom_blank(data=col, aes_string("x","y"))
  else p <- .ggscatter(data = col, x = 'x', y = 'y', 
                       col=col.col,  alpha = alpha.col, 
                       alpha.limits = alpha.limits, shape = shape.col, 
                       geom = geom, lab = lab$col, labelsize = labelsize,
                       pointsize = pointsize, jitter = jitter)
  
  # Add supplementary cols
  if(inherits(X, c('CA', 'ca')) & !hide$col.sup ){
    col_sup <- .get_supp(X, element = "col.sup", axes = axes,
                         select = select.col)
    if(!is.null(col_sup)){
      colnames(col_sup)[2:3] <-  c("x", "y")
      col_sup <- .scale_ca(col_sup, res.ca = X,  element = "col.sup", 
                           type = map, axes = axes)
    }   

    if(!is.null(col_sup)){
      p <- fviz_add(p, df = col_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.col.sup, shape = shape.col,
                    labelsize = labelsize, addlabel = (lab$col.sup & "text" %in% geom),
                    pointsize = pointsize, jitter = jitter)
    }   
  } 
  
  p <- .fviz_finish(p, X, axes) 
  p
  
}




#' @rdname fviz_ca
#' @export 
fviz_ca_biplot <-function(X,  axes = c(1,2), shape.row = 19, shape.col = 17, 
                       geom=c("point", "text"),
                       label = "all", invisible="none", labelsize=4, pointsize =2,
                       col.col ="red", col.col.sup="darkred",  alpha.col = 1,
                       col.row ="blue", col.row.sup="darkblue",  alpha.row = 1,
                       select.col = list(name = NULL, cos2 = NULL, contrib = NULL),
                       select.row = list(name = NULL, cos2 = NULL, contrib = NULL),
                       map ="symmetric", arrows = c(FALSE, FALSE), 
                       jitter = list(what = "label", width = NULL, height = NULL), ...)
{
  
  if(length(intersect(geom, c("point", "text", "arrow"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  if(length(axes) > 2) stop("axes should be of length 2")
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  col <- facto_summarize(X, element = "col", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(col)[2:3] <-  c("x", "y")
  
  # scale row coords according to the type of map
  col <- .scale_ca(col, res.ca = X,  element = "col", 
                   type = map, axes = axes)
  
  col.all <- col
  if(!is.null(select.col)) col <- .select(col, select.col)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.col %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(col.all[, alpha.col])
  
  geom2 <- geom
  if(arrows[1]==TRUE) geom2 <- setdiff(unique(c(geom2, "arrow")), "point")
  p <- fviz_ca_row(X,  axes = axes, shape.row = shape.row, 
        geom=geom2,
        label = label, invisible = invisible, labelsize=labelsize, pointsize = pointsize,
        col.row =col.row, col.row.sup=col.row.sup,  alpha.row = alpha.row, select.row = select.row,
        map = map, jitter = jitter)
  
  
  # geom for columns
  geom2 <- geom
  if(arrows[2]==TRUE) geom2 <- setdiff(unique(c(geom2, "arrow")), "point")
  
  if(!hide$col){
    p <- .ggscatter(p = p, data = col, x = 'x', y = 'y', 
                    col=col.col,  alpha = alpha.col, 
                    alpha.limits = alpha.limits, shape = shape.col, 
                    geom = geom2, lab = lab$col, labelsize = labelsize,
                    pointsize = pointsize, jitter = jitter)
  }
    
  # Add supplementary cols
  if(inherits(X, c('CA', 'ca')) & !hide$col.sup ){
      col_sup <- .get_supp(X, element = "col.sup", axes = axes,
                         select = select.col)
      if(!is.null(col_sup)){
        colnames(col_sup)[2:3] <-  c("x", "y")
        col_sup <- .scale_ca(col_sup, res.ca = X,  element = "col.sup", 
                             type = map, axes = axes)
      }
    if(!is.null(col_sup)){
      p <- fviz_add(p, df = col_sup[, 2:3, drop = FALSE], geom = geom2,
                    color = col.col.sup, shape = shape.col,
                    labelsize = labelsize, addlabel = (lab$col.sup & "text" %in% geom),
                    pointsize = pointsize, jitter = jitter)
    }   
  }    
      
  
  p + labs(title="CA factor map - Biplot")
  
}

#' @rdname fviz_ca
#' @export
fviz_ca <- function(X, ...){
  fviz_ca_biplot(X, ...)
}

#' @include facto_summarize.R
NULL
#' Visualize the contributions of row/column elements
#' 
#' @description
#' This function can be used to visualize the quality of representation (cos2) of rows/columns 
#' from the results of Principal Component Analysis (PCA), 
#' Correspondence Analysis (CA) and 
#' Multiple Correspondence Analysis (MCA) functions.
#' @param ... not used
#' @inheritParams fviz_cos2
#' @details
#' The function fviz_contrib() creates a barplot of row/column contributions. 
#' A reference dashed line is also shown on the barplot. This reference line 
#' corresponds to the expected value if the contribution where uniform.\cr\cr
#' For a given dimension, any row/column with a contribution above the reference line could be 
#' considered as important in contributing to the dimension.
#' 
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal component analysis
#' # ++++++++++++++++++++++++++
#' data(decathlon2)
#' decathlon2.active <- decathlon2[1:23, 1:10]
#' res.pca <- prcomp(decathlon2.active,  scale = TRUE)
#' 
#' # variable contributions on axis 1
#' fviz_contrib(res.pca, choice="var", axes = 1 )
#' # sorting
#' fviz_contrib(res.pca, choice="var", axes = 1, 
#'            sort.val ="asc")
#'            
#' # select the top 7 contributing variables
#' fviz_contrib(res.pca, choice="var", axes = 1, top = 7 )
#' 
#' # Change theme and color
#' fviz_contrib(res.pca, choice="var", axes = 1,
#'          fill = "lightgray", color = "black") +
#'          theme_minimal() + 
#'          theme(axis.text.x = element_text(angle=45))
#'          
#' # Variable contributions on axis 2
#' fviz_contrib(res.pca, choice="var", axes = 2)
#' # Variable contributions on axes 1 + 2
#' fviz_contrib(res.pca, choice="var", axes = 1:2)
#' 
#' # Contributions of individuals on axis 1
#' fviz_contrib(res.pca, choice="ind", axes = 1)
#' 
#' # Correspondence Analysis
#' # ++++++++++++++++++++++++++
#' # Install and load FactoMineR to compute CA
#' # install.packages("FactoMineR")
#' library("FactoMineR")
#' data("housetasks")
#' res.ca <- CA(housetasks, graph = FALSE)
#' 
#' # Visualize row contributions on axes 1
#' fviz_contrib(res.ca, choice ="row", axes = 1)
#' # Visualize row contributions on axes 1 + 2
#' fviz_contrib(res.ca, choice ="row", axes = 1:2)
#' # Visualize column contributions on axes 1
#' fviz_contrib(res.ca, choice ="col", axes = 1)
#' 
#' # Multiple Correspondence Analysis
#' # +++++++++++++++++++++++++++++++++
#' library(FactoMineR)
#' data(poison)
#' res.mca <- MCA(poison, quanti.sup = 1:2, 
#'               quali.sup = 3:4, graph=FALSE)
#'               
#' # Visualize individual contributions on axes 1
#' fviz_contrib(res.mca, choice ="ind", axes = 1)
#' # Select the top 20
#' fviz_contrib(res.mca, choice ="ind", axes = 1, top = 20)
#' # Visualize variable categorie contributions on axes 1
#' fviz_contrib(res.mca, choice ="var", axes = 1)
#'  
#'  }
#'  @export 
fviz_contrib <- function(X, choice = c("row", "col", "var", "ind"), axes=1,
                   fill="steelblue", color = "steelblue",  
                   sort.val = c("desc", "asc", "none"), top = Inf)
{

  title <- .build_title(choice[1], "Contribution", axes)

  dd <- facto_summarize(X, element = choice, result = "contrib", axes = axes)
  contrib <- dd$contrib
  names(contrib) <-rownames(dd)
  
  # expected Average contribution 
  theo_contrib <- 100/length(contrib)
  if(length(axes) > 1) {
    # Adjust variable contributions by the Dimension eigenvalues
    eig <- get_eigenvalue(X)[axes,1]
    theo_contrib <- sum(theo_contrib*eig)
  }
  
  p <- .ggbarplot(contrib, fill =fill, color = color,
                  sort.value = sort.val[1], top = top,
                  title = title, ylab ="Contributions (%)")+
    geom_hline(yintercept=theo_contrib, linetype=2, color="red")
  
  p 
}


#' @describeIn fviz_contrib deprecated function. Use fviz_contrib()
#' @param sortcontrib see the argument sort.val
#' @export 
fviz_pca_contrib <- function(X, choice = c("var", "ind"), axes=1,
                             fill="steelblue", color = "steelblue",  
                             sortcontrib = c("desc", "asc", "none"), top = Inf,...)
{
  
  warning("The function fviz_pca_contrib() is deprecated. ", 
          "Please use the function fviz_contrib() which can handle outputs ",
          " of PCA, CA and MCA functions.")
  
  p <- fviz_contrib(X = X, choice = choice, axes = axes,
               fill = fill, color = color, sort.val = sortcontrib,
               top = top)
  p 
}

#' @include facto_summarize.R
NULL
#' Visualize the quality of representation of rows/columns
#' 
#' @description
#' This function can be used to visualize the quality of representation (cos2) of rows/columns 
#' from the results of Principal Component Analysis (PCA), 
#' Correspondence Analysis (CA) and 
#' Multiple Correspondence Analysis (MCA) functions.
#' @param X an object of class PCA, CA and MCA [FactoMineR]; prcomp and princomp [stats]; 
#'  dudi, pca, coa and acm [ade4]; ca [ca package].
#' @param choice allowed values are "row" and "col" for CA;  "var" and "ind" for PCA or MCA
#' @param axes a numeric vector specifying the dimension(s) of interest.
#' @param fill a fill color for the bar plot.
#' @param color an outline color for the bar plot.
#' @param sort.val a string specifying whether the value should be sorted. 
#' Allowed values are "none" (no sorting), "asc" (for ascending) or "desc" (for descending).
#' @param top a numeric value specifying the number of top elements to be shown.
#'  
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal component analysis
#' # ++++++++++++++++++++++++++
#' data(decathlon2)
#' decathlon2.active <- decathlon2[1:23, 1:10]
#' res.pca <- prcomp(decathlon2.active,  scale = TRUE)
#' 
#' # variable cos2 on axis 1
#' fviz_cos2(res.pca, choice="var", axes = 1 )
#' # sorting
#' fviz_cos2(res.pca, choice="var", axes = 1, 
#'            sort.val ="asc")
#'            
#' # select the top 7 contributing variables
#' fviz_cos2(res.pca, choice="var", axes = 1, top = 7 )
#' 
#' # Change theme and color
#' fviz_cos2(res.pca, choice="var", axes = 1,
#'          fill = "lightgray", color = "black") +
#'          theme_minimal() + 
#'          theme(axis.text.x = element_text(angle=45))
#'          
#' # Variable cos2 on axis 2
#' fviz_cos2(res.pca, choice="var", axes = 2)
#' # Variable cos2 on axes 1 + 2
#' fviz_cos2(res.pca, choice="var", axes = 1:2)
#' 
#' # cos2 of individuals on axis 1
#' fviz_cos2(res.pca, choice="ind", axes = 1)
#' 
#' # Correspondence Analysis
#' # ++++++++++++++++++++++++++
#' # Install and load FactoMineR to compute CA
#' # install.packages("FactoMineR")
#' library("FactoMineR")
#' data("housetasks")
#' res.ca <- CA(housetasks, graph = FALSE)
#' 
#' # Visualize row cos2 on axes 1
#' fviz_cos2(res.ca, choice ="row", axes = 1)
#' # Visualize row cos2 on axes 1 + 2
#' fviz_cos2(res.ca, choice ="row", axes = 1:2)
#' # Visualize column cos2 on axes 1
#' fviz_cos2(res.ca, choice ="col", axes = 1)
#' 
#' # Multiple Correspondence Analysis
#' # +++++++++++++++++++++++++++++++++
#' library(FactoMineR)
#' data(poison)
#' res.mca <- MCA(poison, quanti.sup = 1:2, 
#'               quali.sup = 3:4, graph=FALSE)
#'               
#' # Visualize individual cos2 on axes 1
#' fviz_cos2(res.mca, choice ="ind", axes = 1)
#' # Select the top 20
#' fviz_cos2(res.mca, choice ="ind", axes = 1, top = 20)
#' # Visualize variable categorie cos2 on axes 1
#' fviz_cos2(res.mca, choice ="var", axes = 1)
#'  }
#'  @export 
fviz_cos2 <- function(X, choice = c("row", "col", "var", "ind"), axes=1,
                   fill="steelblue", color = "steelblue",  
                   sort.val = c("desc", "asc", "none"), top = Inf)
{

   title <- .build_title(choice[1], "Cos2", axes)
  
   dd <- facto_summarize(X, element = choice, result = "cos2", axes = axes)
   cos2 <- dd$cos2
   names(cos2) <-rownames(dd)
  p <- .ggbarplot(cos2, fill =fill, color = color,
                  sort.value = sort.val, top = top,
                  title = title, ylab ="Cos2 - Quality of representation")
  
  p 
}

#' @include get_mca.R
 NULL
#' Visualize Multiple Correspondence Analysis
#' 
#' @description
#' Graph of individuals/variables from the output of Multiple Correspondence Analysis (MCA).\cr\cr
#' \itemize{
#' \item{fviz_mca_ind(): Graph of individuals}
#' \item{fviz_mca_var(): Graph of variables}
#' \item{fviz_mca_biplot(): Biplot of individuals and variables}
#' \item{fviz_mca(): An alias of fviz_mca_biplot()}
#' }
#' @param X an object of class MCA [FactoMineR], acm [ade4].
#' @inheritParams fviz_pca
#' @param label a text specifying the elements to be labelled.
#'  Default value is "all".
#'  Allowed values are "none" or the combination of c("ind", "ind.sup","var", "quali.sup",  "quanti.sup"). 
#'  "ind" can be used to label only active individuals. 
#'  "ind.sup" is for supplementary individuals.
#' "var" is for active variable categories.
#'  "quali.sup" is for supplementary qualitative variable categories. 
#' "quanti.sup" is for quantitative supplementary variables.
#' @param invisible a text specifying the elements to be hidden on the plot.
#'  Default value is "none".
#'  Allowed values are the combination of c("ind", "ind.sup","var", "quali.sup",  "quanti.sup").
#' @param habillage an optional factor variable for coloring
#'  the observations by groups. Default value is "none".
#'  If X is an MCA object from FactoMineR package, habillage can also specify
#'  the index of the factor variable in the data. 
#' @param col.ind,col.var color for individuals and variables, respectively.
#'  Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the colors for individuals/variables are automatically controlled by their qualities ("cos2"),
#'  contributions ("contrib"), coordinates (x^2 + y^2 , "coord"), x values("x") or y values("y").
#'  To use automatic coloring (by cos2, contrib, ....), make sure that habillage ="none".
#' @param alpha.ind,alpha.var controls the transparency of
#'  individual and variable colors, respectively.
#' The value can variate from 0 (total transparency) to 1 (no transparency).
#' Default value is 1. Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the transparency for individual/variable colors are automatically controlled by their qualities ("cos2"),
#'  contributions ("contrib"), coordinates (x^2 + y^2 , "coord"), x values("x") or y values("y").
#'  To use this, make sure that habillage ="none".
#' @param shape.ind,shape.var point shapes of individuals and variables
#' @param col.quanti.sup,col.quali.sup a color for the quantitative/qualitative supplementary variables.
#' @param select.ind,select.var a selection of individuals/variables to be drawn. 
#' Allowed values are NULL or a list containing the arguments name, cos2 or contrib: 
#' \itemize{
#' \item name is a character vector containing individuals/variables to be drawn
#' \item cos2 if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn. 
#' if cos2 > 1, ex: 5, then the top 5 individuals/variables with the highest cos2 are drawn.
#' \item contrib if contrib > 1, ex: 5,  then the top 5 individuals/variables with the highest cos2 are drawn
#' }
#' @param ... Arguments to be passed to the function fviz_mca_biplot()
#' @param map character string specifying the map type. Allowed options include: 
#' "symmetric", "rowprincipal", "colprincipal", "symbiplot", "rowgab", 
#' "colgab", "rowgreen" and "colgreen". See details
#' @param arrows Vector of two logicals specifying if the plot should contain
#'  points (FALSE, default) or arrows (TRUE).
#'  First value sets the rows and the second value sets the columns.
#' @param jitter a parameter used to jitter the points in order to reduce overplotting. 
#' It's a list containing the objects what, width and height (i.e jitter = list(what, width, height)). 
#' \itemize{
#' \item what: the element to be jittered. Possible values are "point" or "p"; "label" or "l"; "both" or "b"
#' \item width: degree of jitter in x direction
#' \item height: degree of jitter in y direction
#' }
#' @details The default plot of MCA is a "symmetric" plot in which both rows and 
#' columns are in principal coordinates. In this situation, it's not possible 
#' to interpret the distance between row points and column points. To overcome this 
#' problem, the simplest way is to make an asymmetric plot. This means that, 
#' the column profiles must be presented in row space or vice-versa. 
#' The allowed options for the argument map are: 
#' \itemize{
#' \item "rowprincipal" or "colprincipal": asymmetric plots with either rows in principal 
#' coordinates and columns in standard coordinates, or vice versa. 
#' These plots preserve row metric or column metric respectively.
#' \item "symbiplot": Both rows and columns are scaled to have variances 
#' equal to the singular values (square roots of eigenvalues), 
#' which gives a symmetric biplot but does not preserve row or column metrics. 
#' \item "rowgab" or "colgab": Asymmetric maps, proposed by Gabriel & Odoroff (1990), 
#' with rows (respectively, columns) in 
#' principal coordinates and columns (respectively, rows) in standard coordinates 
#' multiplied by the mass of the corresponding point.
#' \item "rowgreen" or "colgreen": The so-called contribution biplots 
#' showing visually the most contributing points (Greenacre 2006b). 
#' These are similar to "rowgab" and "colgab" except that the points 
#' in standard coordinates are multiplied by the square root of the corresponding masses, 
#' giving reconstructions of the standardized residuals.
#' }
#'  
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Multiple Correspondence Analysis
#' # ++++++++++++++++++++++++++++++
#' # Install and load FactoMineR to compute MCA
#' # install.packages("FactoMineR")
#' library("FactoMineR")
#' data(poison)
#' poison.active <- poison[1:55, 5:15]
#' head(poison.active)
#' res.mca <- MCA(poison.active, graph=FALSE)
#' 
#' # Graph of individuals
#' # +++++++++++++++++++++
#' # Default plot
#' fviz_mca_ind(res.mca)
#' # Change title and axis labels
#' fviz_mca_ind(res.mca) +
#'  labs(title = "MCA", x = "Dim.1", y ="Dim.2" )
#' # Change axis limits by specifying the min and max
#' fviz_mca_ind(res.mca) + 
#'    xlim(-0.8, 1.5) + ylim (-1.5, 1.5)
#' # Use text only
#' fviz_mca_ind(res.mca, geom = "text")
#' # Use points only
#' fviz_mca_ind(res.mca, geom="point")
#' # Change the size of points
#' fviz_mca_ind(res.mca, geom="point", pointsize = 4)
#' # Change point color and theme
#' fviz_mca_ind(res.mca, col.ind = "blue")+
#'    theme_minimal()
#' # Reduce overplotting
#' fviz_mca_ind(res.mca, jitter = list(width = 0.2, height = 0.2))
#'    
#' # Control automatically the color of individuals 
#' # using the cos2 or the contributions
#' # cos2 = the quality of the individuals on the factor map
#' fviz_mca_ind(res.mca, col.ind="cos2") 
#' # Gradient color
#' fviz_mca_ind(res.mca, col.ind="cos2") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=0.4)
#' # Change the theme and use only points
#' fviz_mca_ind(res.mca, col.ind="cos2", geom = "point") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=0.4)+ theme_minimal()
#'       
#' # Color by the contributions   
#' fviz_mca_ind(res.mca, col.ind="contrib") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=1.5)
#'       
#' # Control the transparency of the color by the
#' # contributions
#' fviz_mca_ind(res.mca, alpha.ind="contrib") +
#'      theme_minimal()        
#'              
#' # Color individuals by groups
#' grp <- as.factor(poison.active[, "Vomiting"])
#' fviz_mca_ind(res.mca, label="none", habillage=grp)
#' # Add ellipses
#' p <- fviz_mca_ind(res.mca, label="none", habillage=grp, 
#'              addEllipses=TRUE, ellipse.level=0.95)
#' print(p)
#' # Change group colors using RColorBrewer color palettes
#' p + scale_color_brewer(palette="Dark2") +
#'    theme_minimal()
#' p + scale_color_brewer(palette="Paired") +
#'      theme_minimal()
#' p + scale_color_brewer(palette="Set1") +
#'      theme_minimal()
#'              
#' # Select and visualize individuals with cos2 >= 0.4
#' fviz_mca_ind(res.mca, select.ind = list(cos2 = 0.4))
#' # Select the top 20 according to the cos2
#' fviz_mca_ind(res.mca, select.ind = list(cos2 = 20))
#' # Select the top 20 contributing individuals
#' fviz_mca_ind(res.mca, select.ind = list(contrib = 20))
#' # Select by names
#' fviz_mca_ind(res.mca, 
#' select.ind = list(name = c("44", "38", "53",  "39")))
#' 
#'  
#' # Graph of variable categories
#' # ++++++++++++++++++++++++++++
#' # Default plot
#' fviz_mca_var(res.mca)
#' # Change color and theme
#' fviz_mca_var(res.mca, col.var="steelblue")+
#'  theme_minimal()
#'  
#' # Control variable colors using their contributions
#' fviz_mca_var(res.mca, col.var = "contrib")+
#'  scale_color_gradient2(low = "white", mid = "blue", 
#'            high = "red", midpoint = 2) +
#'  theme_minimal()          
#' # Control the transparency of variables using their contributions
#' fviz_mca_var(res.mca, alpha.var = "contrib") +
#'    theme_minimal()
#'    
#' # Select and visualize categories with cos2 >= 0.4
#' fviz_mca_var(res.mca, select.var = list(cos2 = 0.4))
#' # Select the top 10 contributing variable categories
#' fviz_mca_var(res.mca, select.var = list(contrib = 10))
#' # Select by names
#' fviz_mca_var(res.mca, 
#'  select.var= list(name = c("Courg_n", "Fever_y", "Fever_n")))
#'     
#' # biplot
#' # ++++++++++++++++++++++++++
#' fviz_mca_biplot(res.mca)
#' # Keep only the labels for variable categories
#' fviz_mca_biplot(res.mca, label ="var")
#' # Keep only labels for individuals
#' fviz_mca_biplot(res.mca, label ="ind")
#' # Hide variable categories
#' fviz_mca_biplot(res.mca, invisible ="var")
#' # Hide individuals
#' fviz_mca_biplot(res.mca, invisible ="ind")
#'# Control automatically the color of individuals using the cos2
#' fviz_mca_biplot(res.mca, label ="var", col.ind="cos2") +
#'        theme_minimal()
#' # Change the color by groups, add ellipses
#' fviz_mca_biplot(res.mca, label="var", col.var ="blue",
#'    habillage=grp, addEllipses=TRUE, ellipse.level=0.95) + 
#'    theme_minimal() 
#'                
#' # Select the top 30 contributing individuals
#' # And the top 10 variables
#' fviz_mca_biplot(res.mca,  
#'                select.ind = list(contrib = 30),
#'                select.var = list(contrib = 10)) 
#' 
#'  }
#' @name fviz_mca
#' @rdname fviz_mca 
#' @export 
fviz_mca_ind <- function(X,  axes = c(1,2), geom=c("point", "text"),
                         label = "all", invisible="none", 
                         labelsize=4, pointsize = 2,
                         habillage="none", addEllipses=FALSE, ellipse.level = 0.95,
                         col.ind = "blue", col.ind.sup = "darkblue", alpha.ind =1,
                         shape.ind = 19,
                         select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                         map ="symmetric",
                         jitter = list(what = "label", width = NULL, height = NULL), ...)
{
  
  if(length(intersect(geom, c("point", "text", "arrow"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  if(length(axes) > 2) stop("axes should be of length 2")
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  ind <- facto_summarize(X, element = "ind", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(ind)[2:3] <-  c("x", "y")
  
  # scale ind coords according to the type of map
  ind <- .scale_ca(ind, res.ca = X,  element = "ind", 
                   type = map, axes = axes)
  
  # Selection
  ind.all <- ind
  if(!is.null(select.ind)) ind <- .select(ind, select.ind)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.ind %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(ind.all[, alpha.ind])
  
  # No qualitative variable to color individuals
  if(habillage[1]=="none"){ 
    p <- ggplot() 
    if(hide$ind) p <-ggplot()+geom_blank(data=ind, aes_string('x','y'))
    else p <- .ggscatter(data = ind, x = 'x', y = 'y', 
                         col=col.ind,  alpha = alpha.ind, 
                         alpha.limits = alpha.limits, shape = shape.ind, 
                         geom = geom, lab = lab$ind, labelsize = labelsize,
                         pointsize = pointsize, jitter = jitter)
  }
  # qualitative variable is used to color the individuals
  else{
    
    # Plot individuals
    p <- ggplot()
    if(hide$ind & hide$quali) p <-ggplot()+geom_blank(data=ind, aes_string('x','y'))
    
    if(is.factor(habillage)){ 
      if(nrow(ind)!=length(habillage))
        stop("The number of active individuals used in the MCA is different ",
             "from the length of the factor habillage. Please, remove the supplementary ",
             "individuals in the variable habillage.")
      name.quali <- "Groups"
      ind <- cbind.data.frame(Groups = habillage, ind)
      ind[, 1]<-as.factor(ind[,1])
    }
    # X is from FactoMineR outputs
    else if(inherits(X, "MCA")){
      data <- X$call$X
      if (is.numeric(habillage)) name.quali <- colnames(data)[habillage]
      else name.quali <- habillage 
      ind <- cbind.data.frame(data[rownames(ind),name.quali], ind)
      colnames(ind)[1]<-name.quali
      ind[, 1]<-as.factor(ind[,1])
    }
    
    if(!hide$ind) {
      
      label_coord <- ind
      # jittering
      if(jitter$what %in% c("both", "b")){
        label_coord <- ind <- .jitter(ind, jitter)
      }
      else if(jitter$what %in% c("point", "p")){
        ind <- .jitter(ind, jitter)
      }
      else if(jitter$what %in% c("label", "l")){
        label_coord <- .jitter(label_coord, jitter)
      }
      
      if("point" %in% geom) 
        p <- p+geom_point(data = ind, 
                          aes_string('x', 'y', color=name.quali, shape = name.quali),
                          pointsize = pointsize)
      if(lab$ind & "text" %in% geom) 
        p <- p + geom_text(data = label_coord, 
                           aes_string('x', 'y', label = 'name',
                                      color=name.quali, shape = name.quali),  size = labelsize, vjust = -0.7)
    }
    
    if(!hide$quali){   
      coord_quali.sup <- .get_coord_quali(ind$x, ind$y, groups = ind[,1])
      coord_quali.sup <- cbind.data.frame(name = rownames(coord_quali.sup),
                                          coord_quali.sup)
      colnames(coord_quali.sup)[1] <- name.quali
      coord_quali.sup[, 1] <- as.factor(coord_quali.sup[,1])
      
      if("point" %in% geom) 
        p <- p + geom_point(data=coord_quali.sup,
                            aes_string('x', 'y', color=name.quali, shape=name.quali),
                            size=pointsize*2)    
      if(lab$quali & "text" %in% geom)
        p <- p + geom_text(data=coord_quali.sup, 
                           aes_string('x', 'y', color=name.quali),
                           label=rownames(coord_quali.sup), size=labelsize, vjust=-1)
    }
    if(addEllipses){
      ell <- .get_ellipse_by_groups(ind$x, ind$y,
                                    groups = ind[, name.quali], ellipse.level=ellipse.level)
      colnames(ell)<-c(name.quali, "x", "y")
      ell[, 1]<-as.factor(ell[,1])
      p <- p + geom_path(data = ell, aes_string('x', 'y', color = name.quali, group = name.quali))
    }
    
    
  }
  
  # Add supplementary quantitative individuals
  # Available only in FactoMineR
  if(inherits(X, 'MCA') & !hide$ind.sup){
    ind_sup <- .get_supp(X, element = "ind.sup", axes = axes,
                         select = select.ind)
    if(!is.null(ind_sup)) {
      colnames(ind_sup)[2:3] <-  c("x", "y")
      ind_sup <- .scale_ca(ind_sup, res.ca = X,  element = "ind.sup", 
                           type = map, axes = axes)
    }
    if(!is.null(ind_sup)){
      p <- fviz_add(p, df = ind_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.ind.sup, shape = 19, pointsize = pointsize,
                    labelsize = labelsize, addlabel = (lab$ind.sup & "text" %in% geom), jitter = jitter )
    }  
  }
  
  p <- .fviz_finish(p, X, axes) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_vline(xintercept = 0, color = "black", linetype="dashed") +
    labs(title = "Individuals factor map - MCA")
  
  
  p
}


#' @rdname fviz_mca
#' @export 
fviz_mca_var <- function(X, axes=c(1,2), geom=c("point", "text"), label="all",  invisible ="none",
                         labelsize=4, pointsize = 2, col.var="red", alpha.var=1, shape.var = 17, 
                         col.quanti.sup="blue",  col.quali.sup = "darkgreen", 
                         select.var = list(name = NULL, cos2 = NULL, contrib = NULL),
                         map ="symmetric", jitter = list(what = "label", width = NULL, height = NULL))
{
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <-  c("x", "y")
  
  # scale coords according to the type of map
  var <- .scale_ca(var, res.ca = X,  element = "var", 
                   type = map, axes = axes)
  
  
  # Selection
  var.all <- var
  if(!is.null(select.var)) var <- .select(var, select.var)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.var %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(var.all[, alpha.var])
  
  p <- ggplot()
  
  if(!hide$var){
    p <-.ggscatter(p = p, data = var, x = 'x', y = 'y', 
                   col=col.var,  alpha = alpha.var, 
                   alpha.limits = alpha.limits, 
                   geom = geom, shape = shape.var,
                   lab = lab$var, labelsize = labelsize,
                   pointsize = pointsize, jitter = jitter)
  }
  
  
  # Add supplementary qualitative variable categories
  # Available only in FactoMineR
  if(inherits(X, 'MCA') & !hide$quali.sup ){
    quali_sup <- .get_supp(X, element = "quali.sup", axes = axes,
                            select = select.var)
    if(!is.null(quali_sup)){
      colnames(quali_sup)[2:3] <-  c("x", "y")
      quali_sup <- .scale_ca(quali_sup, res.ca = X,  element = "quali.sup", 
                           type = map, axes = axes)
    }
    if(!is.null(quali_sup)){
      p <- fviz_add(p, df = quali_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.quali.sup, shape = shape.var,
                    labelsize = labelsize, addlabel = (lab$quali.sup),
                    pointsize = pointsize, jitter = jitter)
    }  
    
  }
  
  p <- .fviz_finish(p, X, axes) +
    labs(title = "Variable categories- MCA")
  p 
}



#' @rdname fviz_mca
#' @export
fviz_mca_biplot <- function(X,  axes = c(1,2), geom=c("point", "text"),
                  label = "all", invisible="none", labelsize=4, pointsize = 2,
                  habillage="none", addEllipses=FALSE, ellipse.level = 0.95,
                  col.ind = "blue", col.ind.sup = "darkblue", alpha.ind =1,
                  col.var="red", alpha.var=1, col.quanti.sup="blue",
                  col.quali.sup = "darkgreen", 
                  shape.ind = 19, shape.var = 17, 
                  select.var = list(name = NULL, cos2 = NULL, contrib = NULL),
                  select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                  map ="symmetric", arrows = c(FALSE, FALSE), 
                  jitter = list(what = "label", width = NULL, height = NULL), ...)
{
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <-  c("x", "y")
  
  # scale coords according to the type of map
  var <- .scale_ca(var, res.ca = X,  element = "var", 
                   type = map, axes = axes)
  
  # Selection
  var.all <- var
  if(!is.null(select.var)) var <- .select(var, select.var)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.var %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(var.all[, alpha.var])
  
  
  # Individuals
  geom2 <- geom
  if(arrows[1]==TRUE) geom2 <- setdiff(unique(c(geom2, "arrow")), "point")
  p <- fviz_mca_ind(X,  axes = axes, geom = geom2, label = label, invisible=invisible,
          labelsize=labelsize, pointsize = pointsize,
          col.ind = col.ind, col.ind.sup = col.ind.sup, alpha.ind=alpha.ind,
          shape.ind=shape.ind,
          habillage=habillage, addEllipses=addEllipses, ellipse.level=ellipse.level,
          select.ind = select.ind, jitter = jitter)
    
  # geometry for variable
  geom2 <- geom
  if(arrows[2]==TRUE) geom2 <- setdiff(unique(c(geom2, "arrow")), "point")
  
  if(!hide$var){
    p <-.ggscatter(p = p, data = var, x = 'x', y = 'y', 
                   col=col.var,  alpha = alpha.var, 
                   alpha.limits = alpha.limits, 
                   geom =  geom2, shape = shape.var,
                   lab = lab$var, labelsize = labelsize, pointsize = pointsize, jitter = jitter)
  }
  
  # Add supplementary qualitative variable categories
  # Available only in FactoMineR
  if(inherits(X, 'MCA') & !hide$quali.sup ){
    quali_sup <- .get_supp(X, element = "quali.sup", axes = axes,
                           select = select.var)
    if(!is.null(quali_sup)){
      colnames(quali_sup)[2:3] <- c("x", "y")
      quali_sup <- .scale_ca(quali_sup, res.ca = X,  element = "quali.sup", 
                             type = map, axes = axes)
    }
    if(!is.null(quali_sup)){
      p <- fviz_add(p, df = quali_sup[, 2:3, drop = FALSE], geom = geom2,
                    color = col.quali.sup, shape = shape.var,
                    labelsize = labelsize, addlabel = (lab$quali.sup), pointsize = pointsize, jitter = jitter )
    }  
    
  }
  p+labs(title="MCA factor map - Biplot")
}

#' @rdname fviz_mca
#' @export
fviz_mca <- function(X, ...){
  fviz_mca_biplot(X, ...)
}



#+++++++++++++++++++++
# Helper functions
#+++++++++++++++++++++

#+++++++++++
# .get_ellipse() : Compute the concentration ellipse of the points
#+++++++++++
# x : the coordinates of points. It can be a numeric vector, matrix or data.frame
# y : optional y coordinates of points. y is not required when x
# is a matrix or data.frame
# result is a data.frame containing the x and y coordinates of
# the ellipse. Columns are x, y
.get_ellipse <- function(x, y=NULL, ellipse.level = 0.95) {
  if(class(x)%in% c("matrix", "data.frame")){
    y <- x[,2]
    x <- x[,1]
  }
  sigma <- var(cbind(x, y))
  mu <- c(mean(x), mean(y))
  t <- sqrt(qchisq(ellipse.level, df = 2))
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  data.frame(sweep(circle %*% chol(sigma) * t, 2, mu, FUN = '+'))
}

#+++++++++++
# .get_ellipse() : Compute the concentration ellipse of the points by groups
#+++++++++++
# x : the coordinates of points. It can be a numeric vector, matrix or data.frame
# y : optional y coordinates of points. y is not required when x is a matrix or data.frame
# groups  : a factor variable

# result is a data.frame containing the x and y coordinates of
# the ellipse by groups. Columns are : groups, x, y
.get_ellipse_by_groups <-function(x, y=NULL, groups, ellipse.level = 0.95){
  
  if(class(x)%in% c("matrix", "data.frame")){
    y <- x[,2]
    x <- x[,1]
  }
  groups <-as.factor(groups)
  levs <- levels(groups)
  len <- summary(groups) # number of cases per group
  d <- data.frame(x =x, y = y, groups=groups)
  result <- NULL
  for(i in 1:length(levs)){
    res <- .get_ellipse(d[which(groups==levs[i]),, drop=FALSE], ellipse.level=ellipse.level)
    res <- cbind.data.frame(group=rep(levs[i], nrow(res)), res)
    result <- rbind.data.frame(result,res)
  }
  result
}

# Return the coordinates of groups levels
# x : coordinate of individuals on x axis
# y : coordinate of indiviuals on y axis
.get_coord_quali<-function(x, y, groups){
  data.frame(
    x= tapply(x, groups, mean),
    y = tapply(y, groups, mean)
  )
}

#' @include get_pca.R
 NULL
#' Visualize Principal Component Analysis
#' 
#' @description
#' Graph of individuals/variables from the output of Principal Component Analysis (PCA). \cr\cr
#' \itemize{
#' \item{fviz_pca_ind(): Graph of individuals}
#' \item{fviz_pca_var(): Graph of variables}
#' \item{fviz_pca_biplot(): Biplot of individuals and variables}
#' \item{fviz_pca(): An alias of fviz_pca_biplot()}
#' }
#'  
#' @param X an object of class PCA [FactoMineR]; prcomp and princomp [stats];
#'  dudi and pca [ade4].
#' @param axes a numeric vector of length 2 specifying the dimensions to be plotted.
#' @param geom a text specifying the geometry to be used for the graph. 
#' Allowed values are the combination of c("point", "arrow", "text"). 
#' Use "point" (to show only points); 
#' "text" to show only labels; c("point", "text") or c("arrow", "text") to show both types.
#' @param label a text specifying the elements to be labelled.
#'  Default value is "all".
#'  Allowed values are "none" or the combination of c("ind", "ind.sup", "quali", "var", "quanti.sup").
#'  "ind" can be used to label only active individuals.
#'  "ind.sup" is for supplementary individuals.
#'  "quali" is for supplementary qualitative variables. "var" is for active variables.
#'  "quanti.sup" is for quantitative supplementary variables.
#' @param invisible a text specifying the elements to be hidden on the plot.
#'  Default value is "none".
#'  Allowed values are the combination of c("ind", "ind.sup", "quali", "var", "quanti.sup").
#' @param labelsize font size for the labels
#' @param pointsize the size of points
#' @param habillage an optional factor variable for coloring
#'  the observations by groups. Default value is "none".
#'  If X is a PCA object from FactoMineR package, habillage can also specify
#'  the supplementary qualitative variable (by its index or name) to be used for
#'  coloring individuals by groups (see ?PCA in FactoMineR). 
#' @param addEllipses logical value.
#'  If TRUE, draws ellipses around the individuals when habillage != "none".
#' @param ellipse.level the size of the concentration ellipse in normal probability
#' @param col.ind,col.var color for individuals and variables, respectively.
#'  Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the colors for individuals/variables are automatically controlled by their 
#'  qualities of representation ("cos2"),
#'  contributions ("contrib"), coordinates (x^2+y^2, "coord"), x values ("x") or y values ("y").
#'  To use automatic coloring (by cos2, contrib, ....), make sure that habillage ="none".
#' @param col.ind.sup color for supplementary individuals
#' @param alpha.ind,alpha.var controls the transparency of
#'  individual and variable colors, respectively.
#' The value can variate from 0 (total transparency) to 1 (no transparency).
#' Default value is 1. Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the transparency for the individual/variable colors are automatically controlled by their qualities ("cos2"),
#'  contributions ("contrib"), coordinates (x^2+y^2, "coord"), x values("x") or y values("y").
#'  To use this, make sure that habillage ="none".
#' @param col.quanti.sup a color for the quantitative supplementary variables.
#' @param col.circle a color for the correlation circle.
#' @param select.ind,select.var a selection of individuals/variables to be drawn. 
#' Allowed values are NULL or a list containing the arguments name, cos2 or contrib: 
#' \itemize{
#' \item name: is a character vector containing individuals/variables to be drawn
#' \item cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn. 
#' if cos2 > 1, ex: 5, then the top 5 individuals/variables with the highest cos2 are drawn.
#' \item contrib: if contrib > 1, ex: 5,  then the top 5 individuals/variables with the highest cos2 are drawn
#' }
#' @param jitter a parameter used to jitter the points in order to reduce overplotting. 
#' It's a list containing the objects what, width and height (i.e jitter = list(what, width, height)). 
#' \itemize{
#' \item what: the element to be jittered. Possible values are "point" or "p"; "label" or "l"; "both" or "b".
#' \item width: degree of jitter in x direction
#' \item height: degree of jitter in y direction
#' }
#' @param ... Arguments to be passed to the function fviz_pca_biplot().
#'  
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal component analysis
#' # ++++++++++++++++++++++++++++++
#' data(iris)
#' res.pca <- prcomp(iris[, -5],  scale = TRUE)
#' 
#' # Graph of individuals
#' # +++++++++++++++++++++
#' # Default plot
#' fviz_pca_ind(res.pca)
#' # Change title and axis labels
#' fviz_pca_ind(res.pca) +
#'  labs(title = "PCA", x = "PC1", y ="PC2" )
#' # Change axis limits by specifying the min and max
#' fviz_pca_ind(res.pca) + 
#'    xlim(-4, 4) + ylim (-4, 4)
#' # Use text only
#' fviz_pca_ind(res.pca, geom="text")
#' # Use points only
#' fviz_pca_ind(res.pca, geom="point")
#' # Change the size of points
#' fviz_pca_ind(res.pca, geom="point", pointsize = 4)
#' # Change point color and theme
#' fviz_pca_ind(res.pca, col.ind = "blue")+
#'    theme_minimal()
#'    
#' # Control automatically the color of individuals 
#' # using the cos2 or the contributions
#' # cos2 = the quality of the individuals on the factor map
#' fviz_pca_ind(res.pca, col.ind="cos2") 
#' # Gradient color
#' fviz_pca_ind(res.pca, col.ind="cos2") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=0.6)
#' # Change the theme and use only points
#' fviz_pca_ind(res.pca, col.ind="cos2", geom = "point") + 
#'       scale_color_gradient2(low="blue", mid="white", 
#'       high="red", midpoint=0.6)+ theme_minimal()
#'       
#' # Color by the contributions   
#' fviz_pca_ind(res.pca, col.ind="contrib") + 
#'       scale_color_gradient2(low="blue", mid="white", 
#'       high="red", midpoint=4)
#'       
#' # Control the transparency of the color by the
#' # contributions
#' fviz_pca_ind(res.pca, alpha.ind="contrib") +
#'      theme_minimal()        
#'              
#' # Color individuals by groups
#' fviz_pca_ind(res.pca, label="none", habillage=iris$Species)
#' # Add ellipses
#' p <- fviz_pca_ind(res.pca, label="none", habillage=iris$Species, 
#'              addEllipses=TRUE, ellipse.level=0.95)
#' print(p)
#'              
#' # Change group color using RColorBrewer color palettes
#' p + scale_color_brewer(palette="Dark2") +
#'      theme_minimal()
#' p + scale_color_brewer(palette="Paired") +
#'      theme_minimal()
#' p + scale_color_brewer(palette="Set1") +
#'      theme_minimal()
#' # Change color manually
#' p + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))           
#' # Select and visualize individuals with cos2 >= 0.96
#' fviz_pca_ind(res.pca, select.ind = list(cos2 = 0.96))
#' # Select the top 20 according to the cos2
#' fviz_pca_ind(res.pca, select.ind = list(cos2 = 20))
#' # Select the top 20 contributing individuals
#' fviz_pca_ind(res.pca, select.ind = list(contrib = 20))
#' # Select by names
#' fviz_pca_ind(res.pca, 
#' select.ind = list(name = c("23", "42", "119")))
#' 
#'  
#' # Graph of variables
#' # ++++++++++++++++++++++++++++
#' # Default plot
#' fviz_pca_var(res.pca)
#' # Use points and text
#' fviz_pca_var(res.pca, geom = c("point", "text"))
#' # Change color and theme
#' fviz_pca_var(res.pca, col.var="steelblue")+
#'  theme_minimal()
#'  
#' # Control variable colors using their contributions
#' fviz_pca_var(res.pca, col.var="contrib")+
#'  scale_color_gradient2(low="white", mid="blue", 
#'            high="red", midpoint=96) +
#'  theme_minimal()          
#' # Control the transparency of variables using their contributions
#' fviz_pca_var(res.pca, alpha.var="contrib") +
#'    theme_minimal()
#'    
#' # Select and visualize variables with cos2 >= 0.96
#' fviz_pca_var(res.pca, select.var = list(cos2 = 0.96))
#' # Select the top 3 contributing variables
#' fviz_pca_var(res.pca, select.var = list(contrib = 3))
#' # Select by names
#' fviz_pca_var(res.pca, 
#'    select.var= list(name = c("Sepal.Width", "Petal.Length")))
#'     
#' # biplot
#' # ++++++++++++++++++++++++++
#' fviz_pca_biplot(res.pca)
#' # Keep only the labels for variables
#' fviz_pca_biplot(res.pca, label ="var")
#' # Keep only labels for individuals
#' fviz_pca_biplot(res.pca, label ="ind")
#' # Hide variables
#' fviz_pca_biplot(res.pca, invisible ="var")
#' # Hide individuals
#' fviz_pca_biplot(res.pca, invisible ="ind")
#'# Control automatically the color of individuals using the cos2
#' fviz_pca_biplot(res.pca, label ="var", col.ind="cos2") +
#'        theme_minimal()
#' # Change the color by groups, add ellipses
#' fviz_pca_biplot(res.pca, label="var", habillage=iris$Species,
#'                addEllipses=TRUE, ellipse.level=0.95) 
#'                
#' # Select the top 30 contributing individuals
#' fviz_pca_biplot(res.pca, label="var", 
#'                select.ind = list(contrib = 30)) 
#' 
#'  }
#'  
#' @rdname fviz_pca
#' @export
fviz_pca <- function(X, ...){
  fviz_pca_biplot(X, ...)
}


#' @rdname fviz_pca 
#' @export 
fviz_pca_ind <- function(X,  axes = c(1,2), geom=c("point", "text"),
                         label = "all", invisible="none", labelsize=4, pointsize = 2,
                         habillage="none", addEllipses=FALSE, ellipse.level = 0.95,
                         col.ind = "black", col.ind.sup = "blue", alpha.ind =1,
                         select.ind = list(name = NULL, cos2 = NULL, contrib = NULL), 
                         jitter = list(what = "label", width = NULL, height = NULL),...)
{
  
  if(length(intersect(geom, c("point", "text", "arrow"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  if(length(axes) > 2) stop("axes should be of length 2")
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  ind <- facto_summarize(X, element = "ind", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(ind)[2:3] <-  c("x", "y")
  
  # Selection
  ind.all <- ind
  if(!is.null(select.ind)) ind <- .select(ind, select.ind)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.ind %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(ind.all[, alpha.ind])
  
  # No qualitative variable to color individuals
  if(habillage[1]=="none"){ 
    p <- ggplot() 
    if(hide$ind) p <-ggplot()+geom_blank(data=ind, aes_string('x','y'))
    else p <- .ggscatter(data = ind, x = 'x', y = 'y', 
                         col=col.ind,  alpha = alpha.ind, 
                         alpha.limits = alpha.limits, shape = 19, 
                         geom = geom, lab = lab$ind, labelsize = labelsize,
                         pointsize = pointsize, jitter = jitter)
  }
  
  # qualitative variable is used to color the individuals
  else{
    # Plot individuals
    p <- ggplot()
    if(hide$ind & hide$quali) p <-ggplot()+geom_blank(data=ind, aes_string('x','y'))
    
    if(is.factor(habillage)){ 
      if(nrow(ind)!=length(habillage))
        stop("The number of active individuals used in the PCA is different ",
             "from the length of the factor habillage. Please, remove the supplementary ",
             "individuals in the variable habillage.")
      name.quali <- "Groups"
      ind <- cbind.data.frame(Groups = habillage, ind)
      ind[, 1]<-as.factor(ind[,1])
    }
    # X is from FactoMineR outputs
    else if(inherits(X, "PCA")){
      data <- X$call$X
      if (is.numeric(habillage)) name.quali <- colnames(data)[habillage]
      else name.quali <- habillage 
      ind <- cbind.data.frame(data[rownames(ind),name.quali], ind)
      colnames(ind)[1]<-name.quali
      ind[, 1]<-as.factor(ind[,1])
    }
    
    if(!hide$ind) {   
      
      
      label_coord <- ind
      # jittering
      if(jitter$what %in% c("both", "b")){
        label_coord <- ind <- .jitter(ind, jitter)
      }
      else if(jitter$what %in% c("point", "p")){
        ind <- .jitter(ind, jitter)
      }
      else if(jitter$what %in% c("label", "l")){
        label_coord <- .jitter(label_coord, jitter)
      }
      
      if("point" %in% geom) 
        p <- p+geom_point(data = ind, 
                          aes_string('x', 'y', color=name.quali, shape = name.quali),
                          size = pointsize)
      if(lab$ind & "text" %in% geom) 
        p <- p + geom_text(data = label_coord, 
                           aes_string('x', 'y', label = 'name',
                                      color=name.quali, shape = name.quali),  size = labelsize, vjust = -0.7)
    }
    
    if(!hide$quali){   
      coord_quali.sup <- .get_coord_quali(ind$x, ind$y, groups = ind[,1])
      coord_quali.sup <- cbind.data.frame(name = rownames(coord_quali.sup),
                                          coord_quali.sup)
      colnames(coord_quali.sup)[1] <- name.quali
      coord_quali.sup[, 1] <- as.factor(coord_quali.sup[,1])
      
      if("point" %in% geom) 
        p <- p + geom_point(data=coord_quali.sup,
                            aes_string('x', 'y', color=name.quali, shape=name.quali),
                            size=pointsize*2)    
      if(lab$quali & "text" %in% geom)
        p <- p + geom_text(data=coord_quali.sup, 
                           aes_string('x', 'y', color=name.quali),
                           label=rownames(coord_quali.sup), size=labelsize, vjust=-1)
    }
    if(addEllipses){
      ell <- .get_ellipse_by_groups(ind$x, ind$y,
                                    groups = ind[, name.quali], ellipse.level=ellipse.level)
      colnames(ell)<-c(name.quali, "x", "y")
      ell[, 1]<-as.factor(ell[,1])
      p <- p + geom_path(data = ell, aes_string('x', 'y', color = name.quali, group = name.quali))
    }
    
    
  }
  
  # Add supplementary quantitative individuals
  # Available only in FactoMineR
  if(inherits(X, 'PCA') & !hide$ind.sup){
    ind_sup <- .get_supp(X, element = "ind.sup", axes = axes,
                         select = select.ind)
    if(!is.null(ind_sup)) colnames(ind_sup)[2:3] <-  c("x", "y")
    if(!is.null(ind_sup)){
      p <- fviz_add(p, df = ind_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.ind.sup, shape = 19, pointsize = pointsize,
                    labelsize = labelsize, addlabel = (lab$ind.sup & "text" %in% geom) , jitter = jitter)
    }  
  }
  
  p <- .fviz_finish(p, X, axes) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_vline(xintercept = 0, color = "black", linetype="dashed") +
    labs(title = "Individuals factor map - PCA")
  
  
  p
}


#' @rdname fviz_pca
#' @export 
fviz_pca_var <- function(X, axes=c(1,2), geom=c("arrow", "text"), 
                         label="all",  invisible ="none",
                         labelsize=4, col.var="black", alpha.var=1, 
                         col.quanti.sup="blue", col.circle ="grey70",
                         select.var = list(name = NULL, cos2 = NULL, contrib = NULL),
                         jitter = list(what = "label", width = NULL, height = NULL))
{
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  scale.unit <- .get_scale_unit(X)
  
  # data frame to be used for plotting
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <-  c("x", "y")
  
  # Selection
  var.all <- var
  if(!is.null(select.var)) var <- .select(var, select.var)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.var %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(var.all[, alpha.var])
  
  # Draw correlation circle
  if(scale.unit){
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- data.frame(xcircle = cos(theta), ycircle = sin(theta))
    p <- ggplot(data = circle, aes_string("xcircle", "ycircle")) +
      geom_path(color=col.circle)+
      geom_hline(yintercept = 0, linetype="dashed")+
      geom_vline(xintercept = 0, linetype="dashed")    
  }
  else p <- ggplot()
  
  if(!hide$var){
    p <-.ggscatter(p = p, data = var, x = 'x', y = 'y', 
                   col=col.var,  alpha = alpha.var, 
                   alpha.limits = alpha.limits, 
                   geom =  geom,
                   lab = lab$var, labelsize = labelsize, jitter = jitter)
  }
  
  # Add supplementary quantitative variables
  # Available only in FactoMineR
  if(inherits(X, 'PCA') & !hide$quanti ){
    
    quanti_sup <- .get_supp(X, element = "quanti", axes = axes,
                            select = select.var)
    if(!is.null(quanti_sup)) colnames(quanti_sup)[2:3] <-  c("x", "y")
    if(!is.null(quanti_sup)){
      p <- fviz_add(p, df = quanti_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.quanti.sup, linetype = 2,
                    labelsize = labelsize, addlabel = (lab$quanti), jitter = jitter )
    }  
    
  }
  
  p <- .fviz_finish(p, X, axes) +
    labs(title = "Variables factor map - PCA")
  p 
}



#' @rdname fviz_pca
#' @export
fviz_pca_biplot <- function(X,  axes = c(1,2), geom=c("point", "text"),
                  label = "all", invisible="none", labelsize=4, pointsize = 2,
                  habillage="none", addEllipses=FALSE, ellipse.level = 0.95,
                  col.ind = "black", col.ind.sup = "blue", alpha.ind =1,
                  col.var="steelblue",  alpha.var=1, col.quanti.sup="blue",
                  col.circle ="grey70", 
                  select.var = list(name = NULL, cos2 = NULL, contrib = NULL),
                  select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                  jitter = list(what = "label", width = NULL, height = NULL), ...)
{
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  scale.unit <- .get_scale_unit(X)
  
  # data frame to be used for plotting
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <-  c("x", "y")
  
  # Selection
  var.all <- var
  if(!is.null(select.var)) var <- .select(var, select.var)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.var %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(var.all[, alpha.var])
  
  pca.ind <- get_pca_ind(X)
  ind <- data.frame(pca.ind$coord[, axes, drop=FALSE])
  colnames(ind)<- c("x", "y")
  
  # rescale variable coordinates
  r <- min(
    (max(ind[,"x"])-min(ind[,"x"])/(max(var[,"x"])-min(var[,"x"]))),
    (max(ind[,"y"])-min(ind[,"y"])/(max(var[,"y"])-min(var[,"y"])))
  )
  var[, c("x", "y")] <- var[, c("x", "y")]*r*0.7
  
  # Individuals
  p <- fviz_pca_ind(X,  axes = axes, geom = geom, label = label, invisible=invisible,
          labelsize=labelsize, pointsize = pointsize,
          col.ind = col.ind, col.ind.sup = col.ind.sup, alpha.ind=alpha.ind,
          habillage=habillage, addEllipses=addEllipses, ellipse.level=ellipse.level,
          select.ind = select.ind, jitter = jitter)

  if(!hide$var){
    p <-.ggscatter(p = p, data = var, x = 'x', y = 'y', 
                   col=col.var,  alpha = alpha.var, 
                   alpha.limits = alpha.limits, 
                   geom =  c("arrow", "text"),
                   lab = lab$var, labelsize = labelsize, jitter = jitter)
  }
  
  # Add supplementary quantitative variables
  # Available only in FactoMineR
  if(inherits(X, 'PCA') & !hide$quanti ){
    quanti_sup <- .get_supp(X, element = "quanti", axes = axes,
                            select = select.var)
    if(!is.null(quanti_sup)) colnames(quanti_sup)[2:3] <-  c("x", "y")
    if(!is.null(quanti_sup)){
      p <- fviz_add(p, df = quanti_sup[, 2:3, drop = FALSE]*r*0.7, geom = c("arrow", "text"),
                    color = col.quanti.sup, linetype = 2,
                    labelsize = labelsize, addlabel = (lab$quanti), jitter = jitter )
    }  
  }
  p+labs(title="Biplot of variables and individuals")
}




#+++++++++++++++++++++
# Helper functions
#+++++++++++++++++++++

#+++++++++++
# .get_ellipse() : Compute the concentration ellipse of the points
#+++++++++++
# x : the coordinates of points. It can be a numeric vector, matrix or data.frame
# y : optional y coordinates of points. y is not required when x
# is a matrix or data.frame
# result is a data.frame containing the x and y coordinates of
# the ellipse. Columns are x, y
.get_ellipse <- function(x, y=NULL, ellipse.level = 0.95) {
  if(class(x)%in% c("matrix", "data.frame")){
    y <- x[,2]
    x <- x[,1]
  }
  sigma <- var(cbind(x, y))
  mu <- c(mean(x), mean(y))
  t <- sqrt(qchisq(ellipse.level, df = 2))
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  data.frame(sweep(circle %*% chol(sigma) * t, 2, mu, FUN = '+'))
}

#+++++++++++
# .get_ellipse() : Compute the concentration ellipse of the points by groups
#+++++++++++
# x : the coordinates of points. It can be a numeric vector, matrix or data.frame
# y : optional y coordinates of points. y is not required when x is a matrix or data.frame
# groups  : a factor variable

# result is a data.frame containing the x and y coordinates of
# the ellipse by groups. Columns are : groups, x, y
.get_ellipse_by_groups <-function(x, y=NULL, groups, ellipse.level = 0.95){
  
  if(class(x)%in% c("matrix", "data.frame")){
    y <- x[,2]
    x <- x[,1]
  }
  groups <-as.factor(groups)
  levs <- levels(groups)
  len <- summary(groups) # number of cases per group
  d <- data.frame(x =x, y = y, groups=groups)
  result <- NULL
  for(i in 1:length(levs)){
    res <- .get_ellipse(d[which(groups==levs[i]),, drop=FALSE], ellipse.level=ellipse.level)
    res <- cbind.data.frame(group=rep(levs[i], nrow(res)), res)
    result <- rbind.data.frame(result,res)
  }
  result
}

# Return the coordinates of groups levels
# x : coordinate of individuals on x axis
# y : coordinate of indiviuals on y axis
.get_coord_quali<-function(x, y, groups){
  data.frame(
    x= tapply(x, groups, mean),
    y = tapply(y, groups, mean)
  )
}



# X : an object of class PCA, princomp, prcomp, dudi
# Return TRUE if the data are scaled to unit variance
.get_scale_unit <-function(X){
  scale_unit <- FALSE
  if(inherits(X, 'PCA')) scale_unit <- X$call$scale.unit
  else if(inherits(X, "prcomp" )) {
    scale_unit <- X$scale
    if(is.numeric(scale_unit)) scale_unit = TRUE
  }
  else if(inherits(X, "princomp")){
    scale_unit <- X$scale
    if(length(unique(scale_unit))>1) scale_unit <- TRUE
    else scale_unit = FALSE
  }
  else if(inherits(X, 'pca') & inherits(X, 'dudi')){
    scale_unit <- X$norm
    if(length(unique(scale_unit))>1) scale_unit <- TRUE
    else scale_unit = FALSE
  }
  else stop("Error in .get_scale_unit function : can't handle an object of class ",
            class(X))
  
  scale_unit
}

#' @include get_pca.R
 NULL
#' Visualize Principal Component Analysis
#' 
#' @description
#' Graph of individuals/variables from the output of Principal Component Analysis (PCA). \cr\cr
#' \itemize{
#' \item{fviz_pca_ind(): Graph of individuals}
#' \item{fviz_pca_var(): Graph of variables}
#' \item{fviz_pca_biplot(): Biplot of individuals and variables}
#' \item{fviz_pca(): An alias of fviz_pca_biplot()}
#' }
#'  
#' @param X an object of class PCA [FactoMineR]; prcomp and princomp [stats];
#'  dudi and pca [ade4].
#' @param axes a numeric vector of length 2 specifying the dimensions to be plotted.
#' @param geom a text specifying the geometry to be used for the graph. 
#' Allowed values are the combination of c("point", "arrow", "text"). 
#' Use "point" (to show only points); 
#' "text" to show only labels; c("point", "text") or c("arrow", "text") to show both types.
#' @param label a text specifying the elements to be labelled.
#'  Default value is "all".
#'  Allowed values are "none" or the combination of c("ind", "ind.sup", "quali", "var", "quanti.sup").
#'  "ind" can be used to label only active individuals.
#'  "ind.sup" is for supplementary individuals.
#'  "quali" is for supplementary qualitative variables. "var" is for active variables.
#'  "quanti.sup" is for quantitative supplementary variables.
#' @param invisible a text specifying the elements to be hidden on the plot.
#'  Default value is "none".
#'  Allowed values are the combination of c("ind", "ind.sup", "quali", "var", "quanti.sup").
#' @param labelsize font size for the labels
#' @param pointsize the size of points
#' @param habillage an optional factor variable for coloring
#'  the observations by groups. Default value is "none".
#'  If X is a PCA object from FactoMineR package, habillage can also specify
#'  the supplementary qualitative variable (by its index or name) to be used for
#'  coloring individuals by groups (see ?PCA in FactoMineR). 
#' @param addEllipses logical value.
#'  If TRUE, draws ellipses around the individuals when habillage != "none".
#' @param ellipse.level the size of the concentration ellipse in normal probability
#' @param col.ind,col.var color for individuals and variables, respectively.
#'  Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the colors for individuals/variables are automatically controlled by their 
#'  qualities of representation ("cos2"),
#'  contributions ("contrib"), coordinates (x^2+y^2, "coord"), x values ("x") or y values ("y").
#'  To use automatic coloring (by cos2, contrib, ....), make sure that habillage ="none".
#' @param col.ind.sup color for supplementary individuals
#' @param alpha.ind,alpha.var controls the transparency of
#'  individual and variable colors, respectively.
#' The value can variate from 0 (total transparency) to 1 (no transparency).
#' Default value is 1. Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#'  In this case, the transparency for the individual/variable colors are automatically controlled by their qualities ("cos2"),
#'  contributions ("contrib"), coordinates (x^2+y^2, "coord"), x values("x") or y values("y").
#'  To use this, make sure that habillage ="none".
#' @param col.quanti.sup a color for the quantitative supplementary variables.
#' @param col.circle a color for the correlation circle.
#' @param select.ind,select.var a selection of individuals/variables to be drawn. 
#' Allowed values are NULL or a list containing the arguments name, cos2 or contrib: 
#' \itemize{
#' \item name: is a character vector containing individuals/variables to be drawn
#' \item cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn. 
#' if cos2 > 1, ex: 5, then the top 5 individuals/variables with the highest cos2 are drawn.
#' \item contrib: if contrib > 1, ex: 5,  then the top 5 individuals/variables with the highest cos2 are drawn
#' }
#' @param jitter a parameter used to jitter the points in order to reduce overplotting. 
#' It's a list containing the objects what, width and height (i.e jitter = list(what, width, height)). 
#' \itemize{
#' \item what: the element to be jittered. Possible values are "point" or "p"; "label" or "l"; "both" or "b".
#' \item width: degree of jitter in x direction
#' \item height: degree of jitter in y direction
#' }
#' @param ... Arguments to be passed to the function fviz_pca_biplot().
#'  
#' @return a ggplot2 plot
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal component analysis
#' # ++++++++++++++++++++++++++++++
#' data(iris)
#' res.pca <- prcomp(iris[, -5],  scale = TRUE)
#' 
#' # Graph of individuals
#' # +++++++++++++++++++++
#' # Default plot
#' fviz_pca_ind(res.pca)
#' # Change title and axis labels
#' fviz_pca_ind(res.pca) +
#'  labs(title = "PCA", x = "PC1", y ="PC2" )
#' # Change axis limits by specifying the min and max
#' fviz_pca_ind(res.pca) + 
#'    xlim(-4, 4) + ylim (-4, 4)
#' # Use text only
#' fviz_pca_ind(res.pca, geom="text")
#' # Use points only
#' fviz_pca_ind(res.pca, geom="point")
#' # Change the size of points
#' fviz_pca_ind(res.pca, geom="point", pointsize = 4)
#' # Change point color and theme
#' fviz_pca_ind(res.pca, col.ind = "blue")+
#'    theme_minimal()
#'    
#' # Control automatically the color of individuals 
#' # using the cos2 or the contributions
#' # cos2 = the quality of the individuals on the factor map
#' fviz_pca_ind(res.pca, col.ind="cos2") 
#' # Gradient color
#' fviz_pca_ind(res.pca, col.ind="cos2") + 
#'       scale_color_gradient2(low="white", mid="blue", 
#'       high="red", midpoint=0.6)
#' # Change the theme and use only points
#' fviz_pca_ind(res.pca, col.ind="cos2", geom = "point") + 
#'       scale_color_gradient2(low="blue", mid="white", 
#'       high="red", midpoint=0.6)+ theme_minimal()
#'       
#' # Color by the contributions   
#' fviz_pca_ind(res.pca, col.ind="contrib") + 
#'       scale_color_gradient2(low="blue", mid="white", 
#'       high="red", midpoint=4)
#'       
#' # Control the transparency of the color by the
#' # contributions
#' fviz_pca_ind(res.pca, alpha.ind="contrib") +
#'      theme_minimal()        
#'              
#' # Color individuals by groups
#' fviz_pca_ind(res.pca, label="none", habillage=iris$Species)
#' # Add ellipses
#' p <- fviz_pca_ind(res.pca, label="none", habillage=iris$Species, 
#'              addEllipses=TRUE, ellipse.level=0.95)
#' print(p)
#'              
#' # Change group color using RColorBrewer color palettes
#' p + scale_color_brewer(palette="Dark2") +
#'      theme_minimal()
#' p + scale_color_brewer(palette="Paired") +
#'      theme_minimal()
#' p + scale_color_brewer(palette="Set1") +
#'      theme_minimal()
#' # Change color manually
#' p + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))           
#' # Select and visualize individuals with cos2 >= 0.96
#' fviz_pca_ind(res.pca, select.ind = list(cos2 = 0.96))
#' # Select the top 20 according to the cos2
#' fviz_pca_ind(res.pca, select.ind = list(cos2 = 20))
#' # Select the top 20 contributing individuals
#' fviz_pca_ind(res.pca, select.ind = list(contrib = 20))
#' # Select by names
#' fviz_pca_ind(res.pca, 
#' select.ind = list(name = c("23", "42", "119")))
#' 
#'  
#' # Graph of variables
#' # ++++++++++++++++++++++++++++
#' # Default plot
#' fviz_pca_var(res.pca)
#' # Use points and text
#' fviz_pca_var(res.pca, geom = c("point", "text"))
#' # Change color and theme
#' fviz_pca_var(res.pca, col.var="steelblue")+
#'  theme_minimal()
#'  
#' # Control variable colors using their contributions
#' fviz_pca_var(res.pca, col.var="contrib")+
#'  scale_color_gradient2(low="white", mid="blue", 
#'            high="red", midpoint=96) +
#'  theme_minimal()          
#' # Control the transparency of variables using their contributions
#' fviz_pca_var(res.pca, alpha.var="contrib") +
#'    theme_minimal()
#'    
#' # Select and visualize variables with cos2 >= 0.96
#' fviz_pca_var(res.pca, select.var = list(cos2 = 0.96))
#' # Select the top 3 contributing variables
#' fviz_pca_var(res.pca, select.var = list(contrib = 3))
#' # Select by names
#' fviz_pca_var(res.pca, 
#'    select.var= list(name = c("Sepal.Width", "Petal.Length")))
#'     
#' # biplot
#' # ++++++++++++++++++++++++++
#' fviz_pca_biplot(res.pca)
#' # Keep only the labels for variables
#' fviz_pca_biplot(res.pca, label ="var")
#' # Keep only labels for individuals
#' fviz_pca_biplot(res.pca, label ="ind")
#' # Hide variables
#' fviz_pca_biplot(res.pca, invisible ="var")
#' # Hide individuals
#' fviz_pca_biplot(res.pca, invisible ="ind")
#'# Control automatically the color of individuals using the cos2
#' fviz_pca_biplot(res.pca, label ="var", col.ind="cos2") +
#'        theme_minimal()
#' # Change the color by groups, add ellipses
#' fviz_pca_biplot(res.pca, label="var", habillage=iris$Species,
#'                addEllipses=TRUE, ellipse.level=0.95) 
#'                
#' # Select the top 30 contributing individuals
#' fviz_pca_biplot(res.pca, label="var", 
#'                select.ind = list(contrib = 30)) 
#' 
#'  }
#'  
#' @rdname fviz_pca
#' @export
fviz_pca <- function(X, ...){
  fviz_pca_biplot(X, ...)
}


#' @rdname fviz_pca 
#' @export 
fviz_pca_ind <- function(X,  axes = c(1,2), geom=c("point", "text"),
                         label = "all", invisible="none", labelsize=4, pointsize = 2,
                         habillage="none", addEllipses=FALSE, ellipse.level = 0.95,
                         col.ind = "black", col.ind.sup = "blue", alpha.ind =1,
                         select.ind = list(name = NULL, cos2 = NULL, contrib = NULL), 
                         jitter = list(what = "label", width = NULL, height = NULL),...)
{
  
  if(length(intersect(geom, c("point", "text", "arrow"))) == 0)
    stop("The specified value(s) for the argument geom are not allowed ")
  if(length(axes) > 2) stop("axes should be of length 2")
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  # data frame to be used for plotting
  ind <- facto_summarize(X, element = "ind", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(ind)[2:3] <-  c("x", "y")
  
  # Selection
  ind.all <- ind
  if(!is.null(select.ind)) ind <- .select(ind, select.ind)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.ind %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(ind.all[, alpha.ind])
  
  # No qualitative variable to color individuals
  if(habillage[1]=="none"){ 
    p <- ggplot() 
    if(hide$ind) p <-ggplot()+geom_blank(data=ind, aes_string('x','y'))
    else p <- .ggscatter(data = ind, x = 'x', y = 'y', 
                         col=col.ind,  alpha = alpha.ind, 
                         alpha.limits = alpha.limits, shape = 19, 
                         geom = geom, lab = lab$ind, labelsize = labelsize,
                         pointsize = pointsize, jitter = jitter)
  }
  
  # qualitative variable is used to color the individuals
  else{
    # Plot individuals
    p <- ggplot()
    if(hide$ind & hide$quali) p <-ggplot()+geom_blank(data=ind, aes_string('x','y'))
    
    if(is.factor(habillage)){ 
      if(nrow(ind)!=length(habillage))
        stop("The number of active individuals used in the PCA is different ",
             "from the length of the factor habillage. Please, remove the supplementary ",
             "individuals in the variable habillage.")
      name.quali <- "Groups"
      ind <- cbind.data.frame(Groups = habillage, ind)
      ind[, 1]<-as.factor(ind[,1])
    }
    # X is from FactoMineR outputs
    else if(inherits(X, "PCA")){
      data <- X$call$X
      if (is.numeric(habillage)) name.quali <- colnames(data)[habillage]
      else name.quali <- habillage 
      ind <- cbind.data.frame(data[rownames(ind),name.quali], ind)
      colnames(ind)[1]<-name.quali
      ind[, 1]<-as.factor(ind[,1])
    }
    
    if(!hide$ind) {   
      
      
      label_coord <- ind
      # jittering
      if(jitter$what %in% c("both", "b")){
        label_coord <- ind <- .jitter(ind, jitter)
      }
      else if(jitter$what %in% c("point", "p")){
        ind <- .jitter(ind, jitter)
      }
      else if(jitter$what %in% c("label", "l")){
        label_coord <- .jitter(label_coord, jitter)
      }
      
      if("point" %in% geom) 
        p <- p+geom_point(data = ind, 
                          aes_string('x', 'y', color=name.quali, shape = name.quali),
                          size = pointsize)
      if(lab$ind & "text" %in% geom) 
        p <- p + geom_text(data = label_coord, 
                           aes_string('x', 'y', label = 'name',
                                      color=name.quali, shape = name.quali),  size = labelsize, vjust = -0.7)
    }
    
    if(!hide$quali){   
      coord_quali.sup <- .get_coord_quali(ind$x, ind$y, groups = ind[,1])
      coord_quali.sup <- cbind.data.frame(name = rownames(coord_quali.sup),
                                          coord_quali.sup)
      colnames(coord_quali.sup)[1] <- name.quali
      coord_quali.sup[, 1] <- as.factor(coord_quali.sup[,1])
      
      if("point" %in% geom) 
        p <- p + geom_point(data=coord_quali.sup,
                            aes_string('x', 'y', color=name.quali, shape=name.quali),
                            size=pointsize*2)    
      if(lab$quali & "text" %in% geom)
        p <- p + geom_text(data=coord_quali.sup, 
                           aes_string('x', 'y', color=name.quali),
                           label=rownames(coord_quali.sup), size=labelsize, vjust=-1)
    }
    if(addEllipses){
      ell <- .get_ellipse_by_groups(ind$x, ind$y,
                                    groups = ind[, name.quali], ellipse.level=ellipse.level)
      colnames(ell)<-c(name.quali, "x", "y")
      ell[, 1]<-as.factor(ell[,1])
      p <- p + geom_path(data = ell, aes_string('x', 'y', color = name.quali, group = name.quali))
    }
    
    
  }
  
  # Add supplementary quantitative individuals
  # Available only in FactoMineR
  if(inherits(X, 'PCA') & !hide$ind.sup){
    ind_sup <- .get_supp(X, element = "ind.sup", axes = axes,
                         select = select.ind)
    if(!is.null(ind_sup)) colnames(ind_sup)[2:3] <-  c("x", "y")
    if(!is.null(ind_sup)){
      p <- fviz_add(p, df = ind_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.ind.sup, shape = 19, pointsize = pointsize,
                    labelsize = labelsize, addlabel = (lab$ind.sup & "text" %in% geom) , jitter = jitter)
    }  
  }
  
  p <- .fviz_finish(p, X, axes) +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_vline(xintercept = 0, color = "black", linetype="dashed") +
    labs(title = "Individuals factor map - PCA")
  
  
  p
}


#' @rdname fviz_pca
#' @export 
fviz_pca_var <- function(X, axes=c(1,2), geom=c("arrow", "text"), 
                         label="all",  invisible ="none",
                         labelsize=4, col.var="black", alpha.var=1, 
                         col.quanti.sup="blue", col.circle ="grey70",
                         select.var = list(name = NULL, cos2 = NULL, contrib = NULL),
                         jitter = list(what = "label", width = NULL, height = NULL))
{
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  scale.unit <- .get_scale_unit(X)
  
  # data frame to be used for plotting
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <-  c("x", "y")
  
  # Selection
  var.all <- var
  if(!is.null(select.var)) var <- .select(var, select.var)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.var %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(var.all[, alpha.var])
  
  # Draw correlation circle
  if(scale.unit){
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- data.frame(xcircle = cos(theta), ycircle = sin(theta))
    p <- ggplot(data = circle, aes_string("xcircle", "ycircle")) +
      geom_path(color=col.circle)+
      geom_hline(yintercept = 0, linetype="dashed")+
      geom_vline(xintercept = 0, linetype="dashed")    
  }
  else p <- ggplot()
  
  if(!hide$var){
    p <-.ggscatter(p = p, data = var, x = 'x', y = 'y', 
                   col=col.var,  alpha = alpha.var, 
                   alpha.limits = alpha.limits, 
                   geom =  geom,
                   lab = lab$var, labelsize = labelsize, jitter = jitter)
  }
  
  # Add supplementary quantitative variables
  # Available only in FactoMineR
  if(inherits(X, 'PCA') & !hide$quanti ){
    
    quanti_sup <- .get_supp(X, element = "quanti", axes = axes,
                            select = select.var)
    if(!is.null(quanti_sup)) colnames(quanti_sup)[2:3] <-  c("x", "y")
    if(!is.null(quanti_sup)){
      p <- fviz_add(p, df = quanti_sup[, 2:3, drop = FALSE], geom = geom,
                    color = col.quanti.sup, linetype = 2,
                    labelsize = labelsize, addlabel = (lab$quanti), jitter = jitter )
    }  
    
  }
  
  p <- .fviz_finish(p, X, axes) +
    labs(title = "Variables factor map - PCA")
  p 
}



#' @rdname fviz_pca
#' @export
fviz_pca_biplot <- function(X,  axes = c(1,2), geom=c("point", "text"),
                  label = "all", invisible="none", labelsize=4, pointsize = 2,
                  habillage="none", addEllipses=FALSE, ellipse.level = 0.95,
                  col.ind = "black", col.ind.sup = "blue", alpha.ind =1,
                  col.var="steelblue",  alpha.var=1, col.quanti.sup="blue",
                  col.circle ="grey70", 
                  select.var = list(name = NULL, cos2 = NULL, contrib = NULL),
                  select.ind = list(name = NULL, cos2 = NULL, contrib = NULL),
                  jitter = list(what = "label", width = NULL, height = NULL), ...)
{
  
  if(is.null(jitter$what)) jitter$what <- "label"
  
  scale.unit <- .get_scale_unit(X)
  
  # data frame to be used for plotting
  var <- facto_summarize(X, element = "var", 
                         result = c("coord", "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <-  c("x", "y")
  
  # Selection
  var.all <- var
  if(!is.null(select.var)) var <- .select(var, select.var)
  
  # elements to be labelled or hidden
  lab <- .label(label)
  hide <- .hide(invisible)
  
  alpha.limits <- NULL
  if(alpha.var %in% c("cos2","contrib", "coord", "x", "y"))
    alpha.limits = range(var.all[, alpha.var])
  
  pca.ind <- get_pca_ind(X)
  ind <- data.frame(pca.ind$coord[, axes, drop=FALSE])
  colnames(ind)<- c("x", "y")
  
  # rescale variable coordinates
  r <- min(
    (max(ind[,"x"])-min(ind[,"x"])/(max(var[,"x"])-min(var[,"x"]))),
    (max(ind[,"y"])-min(ind[,"y"])/(max(var[,"y"])-min(var[,"y"])))
  )
  var[, c("x", "y")] <- var[, c("x", "y")]*r*0.7
  
  # Individuals
  p <- fviz_pca_ind(X,  axes = axes, geom = geom, label = label, invisible=invisible,
          labelsize=labelsize, pointsize = pointsize,
          col.ind = col.ind, col.ind.sup = col.ind.sup, alpha.ind=alpha.ind,
          habillage=habillage, addEllipses=addEllipses, ellipse.level=ellipse.level,
          select.ind = select.ind, jitter = jitter)

  if(!hide$var){
    p <-.ggscatter(p = p, data = var, x = 'x', y = 'y', 
                   col=col.var,  alpha = alpha.var, 
                   alpha.limits = alpha.limits, 
                   geom =  c("arrow", "text"),
                   lab = lab$var, labelsize = labelsize, jitter = jitter)
  }
  
  # Add supplementary quantitative variables
  # Available only in FactoMineR
  if(inherits(X, 'PCA') & !hide$quanti ){
    quanti_sup <- .get_supp(X, element = "quanti", axes = axes,
                            select = select.var)
    if(!is.null(quanti_sup)) colnames(quanti_sup)[2:3] <-  c("x", "y")
    if(!is.null(quanti_sup)){
      p <- fviz_add(p, df = quanti_sup[, 2:3, drop = FALSE]*r*0.7, geom = c("arrow", "text"),
                    color = col.quanti.sup, linetype = 2,
                    labelsize = labelsize, addlabel = (lab$quanti), jitter = jitter )
    }  
  }
  p+labs(title="Biplot of variables and individuals")
}




#+++++++++++++++++++++
# Helper functions
#+++++++++++++++++++++

#+++++++++++
# .get_ellipse() : Compute the concentration ellipse of the points
#+++++++++++
# x : the coordinates of points. It can be a numeric vector, matrix or data.frame
# y : optional y coordinates of points. y is not required when x
# is a matrix or data.frame
# result is a data.frame containing the x and y coordinates of
# the ellipse. Columns are x, y
.get_ellipse <- function(x, y=NULL, ellipse.level = 0.95) {
  if(class(x)%in% c("matrix", "data.frame")){
    y <- x[,2]
    x <- x[,1]
  }
  sigma <- var(cbind(x, y))
  mu <- c(mean(x), mean(y))
  t <- sqrt(qchisq(ellipse.level, df = 2))
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  data.frame(sweep(circle %*% chol(sigma) * t, 2, mu, FUN = '+'))
}

#+++++++++++
# .get_ellipse() : Compute the concentration ellipse of the points by groups
#+++++++++++
# x : the coordinates of points. It can be a numeric vector, matrix or data.frame
# y : optional y coordinates of points. y is not required when x is a matrix or data.frame
# groups  : a factor variable

# result is a data.frame containing the x and y coordinates of
# the ellipse by groups. Columns are : groups, x, y
.get_ellipse_by_groups <-function(x, y=NULL, groups, ellipse.level = 0.95){
  
  if(class(x)%in% c("matrix", "data.frame")){
    y <- x[,2]
    x <- x[,1]
  }
  groups <-as.factor(groups)
  levs <- levels(groups)
  len <- summary(groups) # number of cases per group
  d <- data.frame(x =x, y = y, groups=groups)
  result <- NULL
  for(i in 1:length(levs)){
    res <- .get_ellipse(d[which(groups==levs[i]),, drop=FALSE], ellipse.level=ellipse.level)
    res <- cbind.data.frame(group=rep(levs[i], nrow(res)), res)
    result <- rbind.data.frame(result,res)
  }
  result
}

# Return the coordinates of groups levels
# x : coordinate of individuals on x axis
# y : coordinate of indiviuals on y axis
.get_coord_quali<-function(x, y, groups){
  data.frame(
    x= tapply(x, groups, mean),
    y = tapply(y, groups, mean)
  )
}



# X : an object of class PCA, princomp, prcomp, dudi
# Return TRUE if the data are scaled to unit variance
.get_scale_unit <-function(X){
  scale_unit <- FALSE
  if(inherits(X, 'PCA')) scale_unit <- X$call$scale.unit
  else if(inherits(X, "prcomp" )) {
    scale_unit <- X$scale
    if(is.numeric(scale_unit)) scale_unit = TRUE
  }
  else if(inherits(X, "princomp")){
    scale_unit <- X$scale
    if(length(unique(scale_unit))>1) scale_unit <- TRUE
    else scale_unit = FALSE
  }
  else if(inherits(X, 'pca') & inherits(X, 'dudi')){
    scale_unit <- X$norm
    if(length(unique(scale_unit))>1) scale_unit <- TRUE
    else scale_unit = FALSE
  }
  else stop("Error in .get_scale_unit function : can't handle an object of class ",
            class(X))
  
  scale_unit
}

#' @include print.factoextra.R
NULL
#' Extract the results for rows/columns - CA
#' 
#' @description
#' Extract all the results (coordinates, squared cosine, contributions and inertia) 
#' for the active row/column variables from Correspondence Analysis (CA) outputs.\cr\cr
#' \itemize{
#' \item get_ca(): Extract the results for rows and columns
#' \item get_ca_row(): Extract the results for rows only
#' \item get_ca_col(): Extract the results for columns only
#' }
#' @param res.ca an object of class CA [FactoMineR], ca [ca], coa [ade4];
#'  correspondence [MASS].
#' @param element the element to subset from the output. Possible values are "row" or "col".
#' @return a list of matrices containing the results for the active rows/columns including : 
#' \item{coord}{coordinates for the rows/columns}
#' \item{cos2}{cos2 for the rows/columns}
#' \item{contrib}{contributions of the rows/columns}
#' \item{inertia}{inertia of the rows/columns}
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Install and load FactoMineR to compute CA
#' # install.packages("FactoMineR")
#'  library("FactoMineR")
#'  data("housetasks")
#'  res.ca <- CA(housetasks, graph = FALSE)
#'  
#' # Result for column variables
#'  col <- get_ca_col(res.ca)
#'  col # print
#'  head(col$coord) # column coordinates
#'  head(col$cos2) # column cos2
#'  head(col$contrib) # column contributions
#'  
#' # Result for row variables
#'  row <- get_ca_row(res.ca)
#'  row # print
#'  head(row$coord) # row coordinates
#'  head(row$cos2) # row cos2
#'  head(row$contrib) # row contributions
#'  
#'  # You can also use the function get_ca()
#'  get_ca(res.ca, "row") # Results for rows
#'  get_ca(res.ca, "col") # Results for columns
#'  }
#' @name get_ca
#' 
#' @rdname get_ca
#' @export 
get_ca <- function(res.ca, element = c("row", "col")){
 elmt <- element[1]
 if(elmt =="row") get_ca_row(res.ca)
 else if(elmt == "col") get_ca_col(res.ca)
 else stop("Allowed values for the argument element are: 'row' or 'col'.")
}


#' @rdname get_ca
#' @export
get_ca_col <- function(res.ca){
  # FactoMineR package
  if(inherits(res.ca, "CA")) cols <- res.ca$col
  
  # ca package
  else if(inherits(res.ca, "ca")){
    # principal coord = standard coord X sqrt(eig)
    coord <- t(apply(res.ca$colcoord, 1, "*", res.ca$sv))
    cos2 <- apply(coord^2, 2, "/", res.ca$coldist^2)
    # col.contrib <- res.ca$colmass * col.coord^2/res.ca$sv^2
    cc <- apply(coord^2, 2, "*", res.ca$colmass)
    contrib <- t(apply(cc, 1, "/", res.ca$sv^2)) *100
    inertia <- res.ca$colinertia
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    # remove supplementary points
    index <- setdiff(1:nrow(res.ca$colcoord), res.ca$colsup)
    cols <- list(coord = coord[index, , drop = FALSE], 
                 contrib = contrib[index, , drop = FALSE],
                 cos2 = cos2[index, , drop = FALSE], inertia = inertia[index]) 
  }
  # Mass package
  else if(inherits(res.ca, "correspondence")){
    # principal coord = standard coord X sqrt(eig)
    coord <- t(apply(res.ca$cscore, 1, "*", res.ca$cor))
    # cos2 = coord^2/d^2
    row.sum <- apply(res.ca$Freq, 1, sum)
    col.sum <- apply(res.ca$Freq, 2, sum)
    n <- sum(res.ca$Freq)
    profile <- t(apply(res.ca$Freq, 1, "/", col.sum))
    average.profile <- row.sum/n
    d2 <- apply(profile, 2, 
                function(row.p, av.p){sum(((row.p - av.p)^2)/av.p)}, 
                average.profile)
    cos2 <- apply(coord^2, 2, "/", d2)
    # contrib <- mass * coord^2/eig
    mass <- col.sum/n
    cc <- apply(coord^2, 2, "*", mass)
    contrib <- t(apply(cc, 1, "/", res.ca$cor^2)) *100
    # inertia = mass * d^2
    inertia <- mass * d2
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    cols <- list(coord = coord, contrib = contrib, cos2 = cos2, inertia = inertia)
  }
  # ade4 package
  else if(inherits(res.ca, "coa") & inherits(res.ca, 'dudi')){
    if (!requireNamespace("ade4", quietly = TRUE)) {
      stop("ade4 package needed for this function to work. Please install it.")
    }
    coord <- res.ca$co
    inertia <- ade4::inertia.dudi(res.ca, row.inertia = FALSE, col.inertia = TRUE)
    cos2 <- abs(inertia$col.rel/10000)[, colnames(coord)]
    contrib <- (inertia$col.abs/100)[, colnames(coord)]
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    cols <- list(coord = coord, contrib = contrib, cos2 = cos2, inertia = NA)
  }
  
  else stop("An object of class : ", class(res.ca), 
            " can't be handled by the function get_ca_col()")
  class(cols)<-c("factoextra", "ca_col")
  return(cols)
}

#' @rdname get_ca
#' @export
get_ca_row <- function(res.ca){
  
  # FactoMineR package
  if(inherits(res.ca, "CA")) row <- res.ca$row
  
  # ca package
  else if(inherits(res.ca, "ca")){
    # principal coord = standard coord X sqrt(eig)
    coord <- t(apply(res.ca$rowcoord, 1, "*", res.ca$sv))
    cos2 <- apply(coord^2, 2, "/", res.ca$rowdist^2)
    # contrib <- res.ca$rowmass * coord^2/res.ca$sv^2
    cc <- apply(coord^2, 2, "*", res.ca$rowmass)
    contrib <- t(apply(cc, 1, "/", res.ca$sv^2)) *100
    inertia <- res.ca$rowinertia
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    # remove supplementary points
    index <- setdiff(1:nrow(res.ca$rowcoord), res.ca$rowsup)
    row <- list(coord = coord[index, , drop = FALSE], 
                contrib = contrib[index, , drop = FALSE],
                cos2 = cos2[index, , drop = FALSE], inertia = inertia[index])  
  }
  # Mass package
  else if(inherits(res.ca, "correspondence")){
    # principal coord = standard coord X sqrt(eig)
    coord <- t(apply(res.ca$rscore, 1, "*", res.ca$cor))
    # cos2 = coord^2/d^2
    row.sum <- apply(res.ca$Freq, 1, sum)
    col.sum <- apply(res.ca$Freq, 2, sum)
    n <- sum(res.ca$Freq)
    profile <- res.ca$Freq/row.sum
    average.profile <- col.sum/n
    d2 <- apply(profile, 1, 
                function(row.p, av.p){sum(((row.p - av.p)^2)/av.p)}, 
                average.profile)
    cos2 <- apply(coord^2, 2, "/", d2)
    # contrib <- mass * coord^2/eig
    mass <- row.sum/n
    cc <- apply(coord^2, 2, "*", mass)
    contrib <- t(apply(cc, 1, "/", res.ca$cor^2)) *100
    # inertia = mass * d^2
    inertia <- mass * d2
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    row <- list(coord = coord, contrib = contrib, cos2 = cos2, inertia = inertia)
  }
  
  # ade4 package
  else if(inherits(res.ca, "coa") & inherits(res.ca, 'dudi')){
    if (!requireNamespace("ade4", quietly = TRUE)) {
      stop("ade4 package needed for this function to work. Please install it.")
    }
    coord <- res.ca$li
    inertia <- ade4::inertia.dudi(res.ca, row.inertia = TRUE, col.inertia = FALSE)
    cos2 <- abs(inertia$row.rel/10000)[, colnames(coord)]
    contrib <- (inertia$row.abs/100)[, colnames(coord)]
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    row <- list(coord = coord, contrib = contrib, cos2 = cos2, inertia = NA)
  }
  else stop("An object of class : ", class(res.ca), 
            " can't be handled by the function get_ca_row()")
  class(row)<-c("factoextra", "ca_row")
  return(row)
}

#' @include print.factoextra.R
NULL
#' Extract the results for individuals/variables - MCA
#' 
#' @description
#' Extract all the results (coordinates, squared cosine and contributions) 
#' for the active individuals/variable categories from Multiple Correspondence Analysis (MCA) outputs.\cr\cr
#' \itemize{
#' \item get_mca(): Extract the results for variables and individuals
#' \item get_mca_ind(): Extract the results for individuals only
#' \item get_mca_var(): Extract the results for variables only
#' }
#' @param res.mca an object of class MCA [FactoMineR], acm [ade4].
#' @param element the element to subset from the output. Possible values are "var" or "ind".
#' @return a list of matrices containing the results for the active 
#' individuals/variable categories including : 
#' \item{coord}{coordinates for the individuals/variable categories}
#' \item{cos2}{cos2 for the individuals/variable categories}
#' \item{contrib}{contributions of the individuals/variable categories}
#' \item{inertia}{inertia of the individuals/variable categories}
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Multiple Correspondence Analysis
#' # ++++++++++++++++++++++++++++++
#' # Install and load FactoMineR to compute MCA
#' # install.packages("FactoMineR")
#' library("FactoMineR")
#' data(poison)
#' poison.active <- poison[1:55, 5:15]
#' head(poison.active[, 1:6])
#' res.mca <- MCA(poison.active, graph=FALSE)
#'  
#'  # Extract the results for variable categories
#'  var <- get_mca_var(res.mca)
#'  print(var)
#'  head(var$coord) # coordinates of variables
#'  head(var$cos2) # cos2 of variables
#'  head(var$contrib) # contributions of variables
#'  
#'  # Extract the results for individuals
#'  ind <- get_mca_ind(res.mca)
#'  print(ind)
#'  head(ind$coord) # coordinates of individuals
#'  head(ind$cos2) # cos2 of individuals
#'  head(ind$contrib) # contributions of individuals
#'  
#'  # You can also use the function get_mca()
#'  get_mca(res.mca, "ind") # Results for individuals
#'  get_mca(res.mca, "var") # Results for variable categories
#'  }
#'  
#' @name get_mca
#' 
#' @rdname get_mca
#' @export 
get_mca <- function(res.mca, element = c("var", "ind")){
 elmt <- element[1]
 if(elmt =="var") get_mca_var(res.mca)
 else if(elmt == "ind") get_mca_ind(res.mca)
 else stop("Allowed values for the argument element are: 'var' or 'ind'.")
}


#' @rdname get_mca
#' @export
get_mca_var <- function(res.mca){
  # FactoMineR package
  if(inherits(res.mca, "MCA")) vars <- res.mca$var
  # ade4 package
  else if(inherits(res.mca, "acm") & inherits(res.mca, 'dudi')){
    if (!requireNamespace("ade4", quietly = TRUE)) {
      stop("ade4 package needed for this function to work. Please install it.")
    }
    coord <- res.mca$co
    inertia <- ade4::inertia.dudi(res.mca, row.inertia = FALSE, col.inertia = TRUE)
    cos2 <- abs(inertia$col.rel/10000)[, colnames(coord)]
    contrib <- (inertia$col.abs/100)[, colnames(coord)]
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    vars <- list(coord = coord, contrib = contrib, cos2 = cos2)
  }
  else stop("An object of class : ", class(res.mca), 
            " can't be handled by the function get_mca_var()")
  class(vars)<-c("factoextra", "mca_var")
  return(vars)
}

#' @rdname get_mca
#' @export
get_mca_ind <- function(res.mca){
  # FactoMineR package
  if(inherits(res.mca, "MCA")) ind <- res.mca$ind
  # ade4 package
  else if(inherits(res.mca, "acm") & inherits(res.mca, 'dudi')){
    if (!requireNamespace("ade4", quietly = TRUE)) {
      stop("ade4 package needed for this function to work. Please install it.")
    }
    coord <- res.mca$li
    inertia <- ade4::inertia.dudi(res.mca, row.inertia = TRUE, col.inertia = FALSE)
    cos2 <- abs(inertia$row.rel/10000)[, colnames(coord)]
    contrib <- (inertia$row.abs/100)[, colnames(coord)]
    colnames(coord) <- colnames(cos2) <- colnames(contrib) <- paste0("Dim.", 1:ncol(coord)) 
    ind <- list(coord = coord, contrib = contrib, cos2 = cos2)
  }
  
  else stop("An object of class : ", class(res.mca), 
            " can't be handled by the function get_mca_ind()")
  class(ind)<-c("factoextra", "mca_ind")
  return(ind)
}

#' @include print.factoextra.R utilities.R
NULL
#' Extract the results for individuals/variables - PCA
#' 
#' @description
#' Extract all the results (coordinates, squared cosine, contributions) for 
#' the active individuals/variables from Principal Component Analysis (PCA) outputs.\cr\cr
#' \itemize{
#' \item get_pca(): Extract the results for variables and individuals
#' \item get_pca_ind(): Extract the results for individuals only
#' \item get_pca_var(): Extract the results for variables only
#' }
#' @param res.pca an object of class PCA [FactoMineR]; 
#' prcomp and princomp [stats]; pca, dudi [adea4].
#' @param element the element to subset from the output. Allowed values are 
#' "var" (for active variables) or "ind" (for active individuals).
#' @param ... not used
#' @return a list of matrices containing all the results for the active individuals/variables including: 
#' \item{coord}{coordinates for the individuals/variables}
#' \item{cos2}{cos2 for the individuals/variables}
#' \item{contrib}{contributions of the individuals/variables}
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#' # Principal Component Analysis
#' # +++++++++++++++++++++++++++++
#'  data(iris)
#'  res.pca <- prcomp(iris[, -5],  scale = TRUE)
#'  # Extract the results for individuals
#'  ind <- get_pca_ind(res.pca)
#'  print(ind)
#'  head(ind$coord) # coordinates of individuals
#'  head(ind$cos2) # cos2 of individuals
#'  head(ind$contrib) # contributions of individuals
#'  
#'  # Extract the results for variables
#'  var <- get_pca_var(res.pca)
#'  print(var)
#'  head(var$coord) # coordinates of variables
#'  head(var$cos2) # cos2 of variables
#'  head(var$contrib) # contributions of variables
#'  
#'  # You can also use the function get_pca()
#'  get_pca(res.pca, "ind") # Results for individuals
#'  get_pca(res.pca, "var") # Results for variable categories
#'  }
#' @name get_pca
#' 
#' @rdname get_pca
#' @export 
get_pca <- function(res.pca, element = c("var", "ind")){
  elmt <- element[1]
  if(elmt =="var") get_pca_var(res.pca)
  else if(elmt == "ind") get_pca_ind(res.pca)
  else stop("Allowed values for the argument element are: 'ind' or 'var'.")
}

#' @rdname get_pca
#' @export
get_pca_ind<-function(res.pca, ...){
  
  # FactoMineR package
  if(inherits(res.pca, 'PCA')) ind <- res.pca$ind
  
  # ade4 package
  else if(inherits(res.pca, 'pca') & inherits(res.pca, 'dudi')){  
    ind.coord <- res.pca$li
    # get the original data
    data <- res.pca$tab
    data <- t(apply(data, 1, function(x){x*res.pca$norm} ))
    data <- t(apply(data, 1, function(x){x+res.pca$cent}))
    ind <- .get_pca_ind_results(ind.coord, data, res.pca$eig,
                                res.pca$cent, res.pca$norm)
  }
  
  # stats package
  else if(inherits(res.pca, 'princomp')){  
    ind.coord <- res.pca$scores
    data <- .prcomp_reconst(res.pca)
    ind <- .get_pca_ind_results(ind.coord, data, res.pca$sdev^2,
                                res.pca$center, res.pca$scale)
    
  }
  else if(inherits(res.pca, 'prcomp')){
    ind.coord <- res.pca$x
    data <- .prcomp_reconst(res.pca)
    ind <- .get_pca_ind_results(ind.coord, data, res.pca$sdev^2,
                                res.pca$center, res.pca$scale)
  }
  else stop("An object of class : ", class(res.pca), 
            " can't be handled by the function get_pca_ind()")
  
  class(ind)<-c("factoextra", "pca_ind")
  
  ind
}


#' @rdname get_pca
#' @export 
get_pca_var<-function(res.pca){
  # FactoMineR package
  if(inherits(res.pca, 'PCA')) var <- res.pca$var
  # ade4 package
  else if(inherits(res.pca, 'pca') & inherits(res.pca, 'dudi')){
    var <- .get_pca_var_results(res.pca$co)
  }
  # stats package
  else if(inherits(res.pca, 'princomp')){   
    # Correlation of variables with the principal component
    var_cor_func <- function(var.loadings, comp.sdev){var.loadings*comp.sdev}
    var.cor <- t(apply(res.pca$loadings, 1, var_cor_func, res.pca$sdev))
    var <- .get_pca_var_results(var.cor)
  }
  else if(inherits(res.pca, 'prcomp')){
    # Correlation of variables with the principal component
    var_cor_func <- function(var.loadings, comp.sdev){var.loadings*comp.sdev}
    var.cor <- t(apply(res.pca$rotation, 1, var_cor_func, res.pca$sdev))
    var <- .get_pca_var_results(var.cor)
  }
  else stop("An object of class : ", class(res.pca), 
            " can't be handled by the function get_pca_var()")
  class(var)<-c("factoextra", "pca_var")
  var
}




# Helper functions
#++++++++++++++++++++

# compute all the results for individuals : coord, cor, cos2, contrib
# ind.coord : coordinates of variables on the principal component
# pca.center, pca.scale : numeric vectors corresponding to the pca
# center and scale respectively
# data : the orignal data used during the pca analysis
# eigenvalues : principal component eigenvalues
.get_pca_ind_results <- function(ind.coord, data, eigenvalues, pca.center, pca.scale ){
  
  eigenvalues <- eigenvalues[1:ncol(ind.coord)]
  
  if(pca.center[1] == FALSE) pca.center <- rep(0, ncol(data))
  if(pca.scale[1] == FALSE) pca.scale <- rep(1, ncol(data))
  
  # Compute the square of the distance between an individual and the
  # center of gravity
  getdistance <- function(ind_row, center, scale){
    return(sum(((ind_row-center)/scale)^2))
  }
  d2 <- apply(data, 1,getdistance, pca.center, pca.scale)
  
  # Compute the cos2
  cos2 <- function(ind.coord, d2){return(ind.coord^2/d2)}
  ind.cos2 <- apply(ind.coord, 2, cos2, d2)
  
  # Individual contributions 
  contrib <- function(ind.coord, eigenvalues, n.ind){
    100*(1/n.ind)*(ind.coord^2/eigenvalues)
  }
  ind.contrib <- t(apply(ind.coord, 1, contrib,  eigenvalues, nrow(ind.coord)))
  
  colnames(ind.coord) <- colnames(ind.cos2) <-
    colnames(ind.contrib) <- paste0("Dim.", 1:ncol(ind.coord)) 
  
  rnames <- rownames(ind.coord)
  if(is.null(rnames)) rnames <- as.character(1:nrow(ind.coord))
  rownames(ind.coord) <- rownames(ind.cos2) <- rownames(ind.contrib) <- rnames
  
  # Individuals coord, cos2 and contrib
  ind = list(coord = ind.coord,  cos2 = ind.cos2, contrib = ind.contrib)
  ind
}

# compute all the results for variables : coord, cor, cos2, contrib
# var.coord : coordinates of variables on the principal component
.get_pca_var_results <- function(var.coord){
  
  var.cor <- var.coord # correlation
  var.cos2 <- var.cor^2 # variable qualities 
  
  # variable contributions (in percent)
  # var.cos2*100/total Cos2 of the component
  comp.cos2 <- apply(var.cos2, 2, sum)
  contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
  var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
  
  colnames(var.coord) <- colnames(var.cor) <- colnames(var.cos2) <-
    colnames(var.contrib) <- paste0("Dim.", 1:ncol(var.coord)) 
  
  # Variable coord, cor, cos2 and contrib
  list(coord = var.coord, cor = var.cor, cos2 = var.cos2, contrib = var.contrib)
}

#' Print method for an object of class factoextra
#' 
#' @description
#' Print method for an object of class factoextra
#' @param x an object of class factoextra
#' @param ... further arguments to be passed to print method
#' @author Alboukadel Kassambara \email{alboukadel.kassambara@@gmail.com}
#' @references http://www.sthda.com
#' @examples
#' \donttest{
#'  data(iris)
#'  res.pca <- princomp(iris[, -5],  cor = TRUE)
#'  ind <- get_pca_ind(res.pca, data = iris[, -5])
#'  print(ind)
#'  }
#'  
#' @export
print.factoextra<-function(x, ...){
  if(!inherits(x, "factoextra"))
    stop("Can't handle data of class ", class(x))
  
  if(inherits(x, "pca_ind")){
    cat("Principal Component Analysis Results for individuals\n",
        "===================================================\n")
    res <- array(data="", dim=c(3,2), dimnames=list(1:3, c("Name", "Description")))
    res[1, ] <- c("$coord", "Coordinates for the individuals")
    res[2, ] <- c("$cos2", "Cos2 for the individuals")
    res[3, ] <- c("$contrib", "contributions of the individuals")
    print(res[1:3,], ...)
  }
  
  else if(inherits(x, "pca_var")){
    cat("Principal Component Analysis Results for variables\n",
        "===================================================\n")
    res <- array(data="", dim=c(4,2), dimnames=list(1:4, c("Name", "Description")))
    res[1, ] <- c("$coord", "Coordinates for the variables")
    res[2, ] <- c("$cor", "Correlations between variables and dimensions")
    res[3, ] <- c("$cos2", "Cos2 for the variables")
    res[4, ] <- c("$contrib", "contributions of the variables")
    print(res[1:4,])
  }
  
  else if(inherits(x, "ca_row")){
    cat("Correspondence Analysis - Results for rows\n",
        "===================================================\n")
    res <- array(data="", dim=c(4,2), dimnames=list(1:4, c("Name", "Description")))
    res[1, ] <- c("$coord", "Coordinates for the rows")
    res[2, ] <- c("$cos2", "Cos2 for the rows")
    res[3, ] <- c("$contrib", "contributions of the rows")
    res[4, ] <- c("$inertia", "Inertia of the rows")
    print(res[1:4,])
  }
  
  else if(inherits(x, "ca_col")){
    cat("Correspondence Analysis - Results for columns\n",
        "===================================================\n")
    res <- array(data="", dim=c(4,2), dimnames=list(1:4, c("Name", "Description")))
    res[1, ] <- c("$coord", "Coordinates for the columns")
    res[2, ] <- c("$cos2", "Cos2 for the columns")
    res[3, ] <- c("$contrib", "contributions of the columns")
    res[4, ] <- c("$inertia", "Inertia of the columns")
    print(res[1:4,])
  }
  else if(inherits(x, "mca_ind")){
    cat("Multiple Correspondence Analysis Results for individuals\n",
        "===================================================\n")
    res <- array(data="", dim=c(3,2), dimnames=list(1:3, c("Name", "Description")))
    res[1, ] <- c("$coord", "Coordinates for the individuals")
    res[2, ] <- c("$cos2", "Cos2 for the individuals")
    res[3, ] <- c("$contrib", "contributions of the individuals")
    print(res[1:3,], ...)
  }
  else if(inherits(x, "mca_var")){
    cat("Multiple Correspondence Analysis Results for variables\n",
        "===================================================\n")
    res <- array(data="", dim=c(3,2), dimnames=list(1:3, c("Name", "Description")))
    res[1, ] <- c("$coord", "Coordinates for categories")
    res[2, ] <- c("$cos2", "Cos2 for categories")
    res[3, ] <- c("$contrib", "contributions of categories")
    print(res[1:3,])
  }
  
  
}

#' @include eigenvalue.R
NULL
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_alpha
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom grid arrow
#' @importFrom grid unit
# Check and get the class of the output of a factor analysis
# ++++++++++++++++++++++++++++
# X: an output of factor analysis (PCA, CA, MCA) 
# from different packages (FactoMineR, ade4, ....) 
.get_facto_class <- function(X){
  if(inherits(X, c('PCA', 'princomp', 'prcomp')))
    facto_class ="PCA"
  else if(inherits(X, 'pca') & inherits(X, 'dudi'))
    facto_class ="PCA"
  else if(inherits(X, c("CA", "ca", "coa", "correspondence"))) facto_class="CA"
  else if(inherits(X, c("MCA", "acm"))) facto_class = "MCA"
  else stop("An object of class : ", class(X), 
            " can't be handled by factoextra")   
}

# Get the result for supplementary points
#++++++++++++++++++++++++++
## X: an output of factor analysis (PCA, CA, MCA) 
## element possible values are "row.sup", "col.sup" (CA);
# quanti, ind.sup (PCA); quali.sup (MCA)
## result the result tobe extracted for the element. Possible values are
#  the combination of c("cos2", "contrib", "coord")
## axes a numeric vector specifying the axes of interest. Default values are 1:2
##  for axes 1 and 2
## select: a selection of variables. See the function .select()
.get_supp <- function(X, element = NULL, axes = 1:2, 
                      result = c("coord", "cos2"), select = NULL){
  
  if(inherits(X, c("CA", "PCA", "MCA"))) {
    exprs <- paste0("X$", element)
    elmt <- eval(parse(text=exprs ))
  }
  else if(inherits(X, "ca")){
    if(element == "col.sup") elmt <- .get_ca_col_sup(X)
    else if(element == "row.sup") elmt <- .get_ca_row_sup(X)
  }
  else stop("An object of class : ", class(X), 
            " can't be handled by the function .get_supp()")
  
  # summarize the result
  res = NULL
  if(!is.null(elmt)){
    # check axes
    if(max(axes) > ncol(elmt$coord))
      stop("The value of the argument axes is incorrect. ",
           "The number of axes in the data is: ", ncol(elmt$coord), 
           ". Please try again with axes between 1 - ", ncol(elmt$coord))
    
    # 1.Extract the coordinates x, y and coord
    if("coord" %in% result){
      dd <- data.frame(elmt$coord[, axes, drop=FALSE])
      coord <- apply(dd^2, 1, sum) # x^2 + y2 + ...
      res = cbind(dd, coord = coord)
    }
    
    # 2. Extract the cos2
    if("cos2" %in% result){
      cos2 <- elmt$cos2[, axes]
      if(length(axes) > 1) cos2 <- apply(cos2, 1, sum, na.rm=TRUE)
      res <- cbind(res, cos2 = cos2)
    }
    
    res <- cbind.data.frame(name =rownames(elmt$coord), res)
    
    # selection of variables
    if(!is.null(select)){
      if(!is.null(select$contrib)) res <- NULL # supp points don't have contrib
      else res <- .select(res, select, check = FALSE)
    }
  }
  res
}


# get supplementary columns from ca output (ca package)
.get_ca_col_sup <- function(res.ca){
  # supplementary points
  index <- res.ca$colsup
  cols <- NULL 
  if(length(index) > 0){
    # principal coord = standard coord X sqrt(eig)
    coord <- t(apply(res.ca$colcoord, 1, "*", res.ca$sv))
    cos2 <- apply(coord^2, 2, "/", res.ca$coldist^2)
    cols <- list(coord = coord[index, , drop = FALSE], 
                cos2 = cos2[index, , drop = FALSE]) 
  }
  cols
}

# get supplementary rows from ca output (ca package)
.get_ca_row_sup <- function(res.ca){
  rows <- NULL
  # supplementary points
  index <- res.ca$rowsup
  if(length(index) > 0){
    # principal coord = standard coord X sqrt(eig)
    coord <- t(apply(res.ca$rowcoord, 1, "*", res.ca$sv))
    cos2 <- apply(coord^2, 2, "/", res.ca$rowdist^2)
    rows <- list(coord = coord[index, , drop = FALSE], 
                                     cos2 = cos2[index, , drop = FALSE]) 
  }
  rows
}

# Scale row coordinates depending on the map type
#++++++++++++++++++++++
## row.res, data: an output of facto_summarize containing the
# name, x, y, coord,... of the row variables
## res.ca an object of class CA, MCA
## type: a character string specifying the map type. 
# Allowed values include "symmetric", "rowprincipal", "colprincipal", "rowgab", "colgab", "rowgreen", "colgreen".
## element: possible values are "row", "col", "row.sup", "col.sup",
# "ind", "var", "ind.sup", "quali.sup"
.scale_ca <- function(data, res.ca, element="row", type="symmetric", axes = 1:2)
{
  res <- NULL
  if(!is.null(data)){
  if(element %in% c("row", "ind")) res <- .scale_ca_row(data, res.ca, type, axes)
  else if(element %in% c("col", "var")) res <- .scale_ca_col(data, res.ca, type, axes)
  else if(element %in% c( "row.sup", "ind.sup")) res <- .scale_ca_rowsupp(data, res.ca, type, axes)
  else if(element %in% c("col.sup", "quali.sup")) res <- .scale_ca_colsupp(data, res.ca, type, axes)
  } 
  res
}

.scale_ca_row <- function(row.res, res.ca, type ="symmetric", axes = 1:2){
  data <- row.res
  res <- data
  
  eig <- get_eigenvalue(res.ca)[axes,1]
  sv <- sqrt(eig)
  mass <- .get_ca_mass(res.ca, "row")[as.vector(data$name)]
  
  # rows/columns in principal coordinates
  if(type %in% c("symmetric", "rowprincipal", "rowgab", "rowgreen")) res <- data
  else{
    # get standard coordinates 
    x <- res$x/sqrt(eig[1])
    y <- res$y/sqrt(eig[2])
    
    # standard coordinate
    if(type %in% c("colprincipal")){
      res$x <- x
      res$y <- y
      res$coord <- res$x^2 + res$y^2
    }
    # standard coordinate X mass
    else if(type %in% c("colgab")){
      res$x <- x*mass
      res$y <- y*mass
      res$coord <- res$x^2 + res$y^2
    }
    
    # standard coordinate X sqrt(mass)
    else if(type %in% c("colgreen")){
      res$x <- x*sqrt(mass)
      res$y <- y*sqrt(mass)
      res$coord <- res$x^2 + res$y^2
    }
    else  if(type %in% c("symbiplot")){
      res$x <- x*sqrt(sv[1])
      res$y <- y*sqrt(sv[2])
      res$coord <- res$x^2 + res$y^2
    }
  }
  
  return(res)
}

# Scale ca column
.scale_ca_col <- function(col.res, res.ca, type ="symmetric", axes = 1:2){
  data <- col.res
  res <- data
  
  eig <- get_eigenvalue(res.ca)[axes,1]
  sv <- sqrt(eig)
  mass <- .get_ca_mass(res.ca, "col")[as.vector(data$name)]
  
  # rows/columns in principal coordinates
  if(type %in% c("symmetric", "colprincipal", "colgab", "colgreen")) res <- data
  else{
    # get standard coordinates 
    x <- res$x/sqrt(eig[1])
    y <- res$y/sqrt(eig[2])
    
    # standard coordinate
    if(type %in% c("rowprincipal")){
      res$x <- x
      res$y <- y
      res$coord <- res$x^2 + res$y^2
    }
    # standard coordinate X mass
    else if(type %in% c("rowgab")){
      res$x <- x*mass
      res$y <- y*mass
      res$coord <- res$x^2 + res$y^2
    }
    
    # standard coordinate X sqrt(mass)
    else if(type %in% c("rowgreen")){
      res$x <- x*sqrt(mass)
      res$y <- y*sqrt(mass)
      res$coord <- res$x^2 + res$y^2
    }
    else  if(type %in% c("symbiplot")){
      res$x <- x*sqrt(sv[1])
      res$y <- y*sqrt(sv[2])
      res$coord <- res$x^2 + res$y^2
    }
  }
  
  return(res)
}

# Scale ca row sup
.scale_ca_rowsupp <- function(rowsup.res, res.ca, type ="symmetric", axes = 1:2){
  data <- rowsup.res
  res <- data
  
  eig <- get_eigenvalue(res.ca)[axes,1]
  sv <- sqrt(eig)
  
  # rows/columns in principal coordinates
  if(type %in% c("symmetric", "rowprincipal", "rowgab", "rowgreen")) res <- data
  else{
    # get standard coordinates 
    x <- res$x/sqrt(eig[1])
    y <- res$y/sqrt(eig[2])
    
    # standard coordinate
    if(type %in% c("colprincipal")){
      res$x <- x
      res$y <- y
      res$coord <- res$x^2 + res$y^2
    }
    # standard coordinate X mass
    else if(type %in% c("colgab")) res <- NULL
    
    # standard coordinate X sqrt(mass)
    else if(type %in% c("colgreen")) res <- NULL
    
    else  if(type %in% c("symbiplot")){
      res$x <- x*sqrt(sv[1])
      res$y <- y*sqrt(sv[2])
      res$coord <- res$x^2 + res$y^2
    }
  }
  
  return(res)
}


.scale_ca_colsupp <- function(colsup.res, res.ca, type ="symmetric", axes = 1:2){
  data <- colsup.res
  res <- data
  
  eig <- get_eigenvalue(res.ca)[axes,1]
  sv <- sqrt(eig)
  
  # rows/columns in principal coordinates
  if(type %in% c("symmetric", "colprincipal", "colgab", "cogreen")) res <- data
  else{
    # get standard coordinates 
    x <- res$x/sqrt(eig[1])
    y <- res$y/sqrt(eig[2])
    
    # standard coordinate
    if(type %in% c("rowprincipal")){
      res$x <- x
      res$y <- y
      res$coord <- res$x^2 + res$y^2
    }
    # standard coordinate X mass
    else if(type %in% c("rowgab")) res <- NULL
    # standard coordinate X sqrt(mass)
    else if(type %in% c("rowgreen")) res <- NULL
    else  if(type %in% c("symbiplot")){
      res$x <- x*sqrt(sv[1])
      res$y <- y*sqrt(sv[2])
      res$coord <- res$x^2 + res$y^2
    }
  }
  
  return(res)
}


# get the mass of an element
#++++++++++++++++++++++++++++
# res.ca : CA object
# element: possible values are "row" or "col"
.get_ca_mass <- function(res.ca, element){
  
  if(inherits(res.ca, "ca")){
    if(element == "row"){
    mass <- res.ca$rowmass
    names(mass) <- res.ca$rownames
    }
    else if(element == "col"){
      mass <- res.ca$colmass
      names(mass) <- res.ca$colnames
    }
  }
  # FactoMiner
  else if(inherits(res.ca, c("CA", "MCA"))){
    if(element %in% c("row", "ind")) mass <- res.ca$call$marge.row
    else if(element %in% c("col", "var")) mass <- res.ca$call$marge.col
    
  }
  # Ade4
  else if(inherits(res.ca, c("coa", "acm"))){
    if(element =="row") mass <- res.ca$lw
    else if(element == "col") mass <- res.ca$cw
  }
  # Mass package
  else if(inherits(res.ca, "correspondence")){
    row.sum <- apply(res.ca$Freq, 1, sum)
    col.sum <- apply(res.ca$Freq, 2, sum)
    n <- sum(res.ca$Freq)
    if(element =="row") mass <- row.sum/n
    else if(element =="col") mass <- col.sum/n
  }
  else stop("An object of class : ", class(res.ca), 
            " can't be handled") 
  return(mass)
}


# Reconstruct data: get the original data
# res.pca is an object of class prcomp or princomp
.prcomp_reconst <- function(res.pca){
  if(inherits(res.pca, "prcomp")){
    centered <- is.numeric(res.pca$center)
    scaled <- is.numeric(res.pca$scale)
    
    if(centered & scaled) 
      t(t(res.pca$x %*% t(res.pca$rotation)) * res.pca$scale + res.pca$center)
    else if(centered)
      t(t(res.pca$x %*% t(res.pca$rotation)) + res.pca$center)
    else if(scaled)
      t(t(res.pca$x %*% t(res.pca$rotation)) * res.pca$scale)
    else t(t(res.pca$x %*% t(res.pca$rotation)))
  }
  
  else if(inherits(res.pca, "princomp")){ 
    if(!is.null(res.pca$scores))
      t(t(res.pca$scores %*% t(res.pca$loadings)) * res.pca$scale + res.pca$center)
    else stop("The object res.pca doesn't have the element scores. ", 
              "Please use the function princomp() with the argument ",
              "scores = TRUE.")
  }
  
  else stop("Can't handle an object of class ", class(res.pca))
}



# Build plot title for contribution, cos2 plots
# used in fviz_cos2, fviz_contrib, etc
#+++++++++++++++++++++++
# element: possible values are row, col, var, ind
# varname: the name of the variable name to be plotted
  # possible values are "Cos2", "Contribution"
# axes: a numeric vector specifying the axes of interest
# example: .build_title("row", "Cos2", 1:2)
.build_title <- function(element, varname,  axes){
  if(varname =="cos2") varname = "Cos2"
  else if(varname =="contrib") varname = "Contribution"
  
  if(element=="row") 
    title <- paste0(varname, " of rows to Dim-", paste(axes, collapse="-"))
  else if(element=="col")
    title <- paste0(varname, " of columns to Dim-", paste(axes, collapse="-"))
  else if(element=="var")
    title <- paste0(varname, " of variables to Dim-", paste(axes, collapse="-"))
  else if(element=="ind")
    title <- paste0(varname, " of individuals to Dim-", paste(axes, collapse="-"))
  
  return(title)
}

# Select rows according to a filter
# ++++++++++++++++++++
# d a data frame containing row names, coordinates, cos2, contrib, etc
# filter: a filter to be applied. Allowed values are NULL or a list containing the arguments either
#  name, cos2 or contrib
# - name: is a character vector containing row names of interest
# - cos2: if cos2 is in [0, 1], ex: 0.6, then rows with a cos2 > 0.6 are extracted.
#   if cos2 > 1, ex: 5, then the top 5 rows with the highest cos2 are extracted
# - contrib: if contrib > 1, ex: 5,  then the top 5 rows with the highest cos2 are extracted
# check: if TRUE, check the data after filtering
.select <- function(d, filter = NULL, check= TRUE){
  
  if(!is.null(filter)){
    
    # Filter by name
    if(!is.null(filter$name)){
      name <- filter$name
      common <- intersect(name, d$name)
      diff <- setdiff(name, d$name)
      #if(check & length(common) == 0) stop("Can't find the specified names")
      # if(check & length(diff)!=0) warning("Can't find the the following name(s): ", diff)
      d <- d[common, , drop = FALSE]
    }
    
    # Filter by cos2
    if(!is.null(filter$cos2) & nrow(d) >= 1){
      # case 1 cos2 is in [0, 1]
      # rows with cos2 > value are selected
      if(0 <= filter$cos2 & filter$cos2 <= 1){
        d <- d[which(d$cos2 >= filter$cos2), , drop = FALSE]
        if(check & nrow(d)==0)
          stop("There are no observations with cos2 >=", filter$cos2, 
               ". Please, change the value of cos2 and try again.")
      }
      # case 2 - cos2 > 1 : the top rows are selected 
      else if(filter$cos2 > 1){
        cos2 <- round(filter$cos2)
        d <- d[order(d$cos2, decreasing = TRUE), , drop = FALSE]
        d <- d[1:min(filter$cos2, nrow(d)),, drop = FALSE]
      }
    }
    
    # Filter by contrib: the top rows are selected 
    if(!is.null(filter$contrib) & nrow(d) >= 1){
      contrib <- round(filter$contrib)
      if(contrib < 1) stop("The value of the argument contrib >", 1)
      d <- d[order(d$contrib, decreasing = TRUE), , drop = FALSE]
      d <- d[1:min(contrib, nrow(d)), , drop = FALSE]
    }
    
  }
  
  return (d)
}


# Add supplementary points

# Make a ggplot2 bar plot
#+++++++++++++++++++++++++
# x: a numeric vector
# title, xlab, ylab : labels for the graph
# fill: a fill color for the bar plot
# color: an outline color for the bar plot
# sort.value: a string specifying whether x should be sorted or not. 
# Allowed values are "none" (no sorting), "asc" (for ascending) or "desc" (for descending)
# top a numeric value specifing the top  elements to be shown
.ggbarplot <- function(x,  title ="", xlab ="", ylab="",
                       fill="steelblue", color = "steelblue",  
                       sort.value = c("desc", "asc", "none"), top = Inf ){
  
 
  # top elements
  if(top!=Inf & top < length(x))
    x <- sort(x, decreasing=TRUE)[1:top]
  # sorting
  if(sort.value[1]=="desc") x <- sort(x, decreasing = TRUE)
  else if(sort.value[1]=="asc") x <- sort(x, decreasing = FALSE)
  # bar names
  if(is.null(names(x))) names(x) <- 1:length(x)
  #data frame for ggplot2
  d <- cbind.data.frame(name = factor(names(x), levels = names(x)), val = x)
  # plot
  p <- ggplot(d, aes_string("name", "val")) + 
    geom_bar(stat="identity", fill=fill, color = color) +
    labs(title = title,  x =xlab, y = ylab)+
    theme(axis.text.x = element_text(angle=45), 
                   axis.title.x = element_blank())
  
  return(p)
}


# Make a scatter plot
#++++++++++++++++++++++++++
## p: a ggplot. If not null, points are added on the existing plot
## data: is a data frame of form: name|x|y|coord|cos2|contrib
## col: point color. The default value is "black".
#  Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#  In this case, the colors for variables are automatically controlled by their qualities ("cos2"),
#  contributions ("contrib"), coordinates (x^2+y2, "coord"), x values("x") or y values("y")
## x, y: contains the name of x and y variables
## alpha: controls the transparency of colors.
# The value can variate from 0 (total transparency) to 1 (no transparency).
# Default value is 1. Possible values include also : "cos2", "contrib", "coord", "x" or "y".
## alpha.limits: the range of alpha values. Used only when alpha is 
# a continuous variable(i.e, cos2, contrib, ...)
## shape: point shape
## pointsize: point size
## lab: if TRUE points are labelled
## geom: a character specifying the geometry to be used for the graph.
#  Allowed values are the combination of c("point", "arrow", "text"). Use "point" (to show only points),
#  "text" to show only labels or c("point", "text") to show both types.
## jitter: to avoid overplotting. Possible values for what are "label", "point", "both"
# it's possible to use also the shortcut "l", "p", "b".
.ggscatter <- function(p = NULL, data, x = 'x', y = 'y', col="black",  alpha = 1,
                       alpha.limits = NULL, shape = 19, pointsize = 2, 
                       geom=c("point", "text"), lab = TRUE, labelsize = 4, 
                       jitter = list(what = "label", width = NULL, height = NULL), ...)
  {
  
  data <- as.data.frame(data)
  label_coord <- data
  
  # jittering
  if(jitter$what %in% c("both", "b")){
    label_coord <- data <- .jitter(data, jitter)
  }
  else if(jitter$what %in% c("point", "p")){
    data <- .jitter(data, jitter)
  }
  else if(jitter$what %in% c("label", "l")){
    label_coord <- .jitter(label_coord, jitter)
  }
  
#   uvar = "xxx" # user variable used for coloring
#   # Color is a factor variable (Groups)
#   if(is.factor(col)){ 
#     if(nrow(data)!=length(col))
#       stop("The number of points and the length of the factor variable for coloring are different.")
#     uvar <- "Groups"
#     data<- cbind.data.frame(Groups = col, data)
#     data[, 1]<-as.factor(data[,1])
#   }
#   # Color is a numeric vector
#   else if(is.numeric(col) & length(col) == nrow(data)){
#     uvar <- "Col.var"
#     data<- cbind.data.frame(Col.var = col, data)
#     data[, 1]<-as.factor(data[,1])
#   }
#   else if(!(col %in% c("cos2","contrib", "coord", "x", "y"))){
#     stop("The argument col is incorrect.")
#   }
  
 
  
  if(is.null(p)) p <- ggplot() 
  # The color and the transparency of variables are automatically controlled by
  # their cos2, contrib,  "x" or "y" coordinates
  if(col %in% c("cos2","contrib", "coord", "x", "y") &
            alpha %in% c("cos2","contrib", "coord", "x", "y"))
  {
    if("point" %in% geom) 
      p <- p + geom_point(data = data, 
                          aes_string(x,y, color=col, alpha=alpha),
                          shape=shape, size = pointsize)
    else if("arrow" %in% geom) 
      p <- p + geom_segment(data = data,
                  aes_string(x = 0, y = 0, xend = x, yend = y, color=col, alpha=alpha),
                  arrow = grid::arrow(length = grid::unit(0.2, 'cm')))
  
    if(lab & "text"%in% geom) 
      p <- p + geom_text(data = label_coord, 
                         aes_string(x,y, label = 'name', 
                                    color=col, alpha=alpha), size = labelsize)
    if(!is.null(alpha.limits)) p <- p + scale_alpha(limits = alpha.limits)
  }
  # Only the color is controlled automatically
  else if(col %in% c("cos2","contrib", "coord", "x", "y")){
    
    if("point" %in% geom) 
      p <- p + geom_point(data = data, aes_string(x,y, color=col),
                          shape=shape, alpha=alpha, size=pointsize)
    else if("arrow" %in% geom) 
      p <- p + geom_segment(data = data,
                            aes_string(x = 0, y = 0, xend = x, yend = y, color=col),
                            arrow = grid::arrow(length = grid::unit(0.2, 'cm')), alpha = alpha)
    
    if(lab & "text" %in% geom) 
      p <- p + geom_text(data = label_coord,
                         aes_string(x, y, color=col),
                         label = data$name,  size = labelsize, alpha=alpha, vjust = -0.7)
  }
  
  # Only the transparency is controlled automatically
  else if(alpha %in% c("cos2","contrib", "coord", "x", "y")){
    
    if("point" %in% geom) 
      p <- p + geom_point(data = data, 
                          aes_string(x, y, alpha=alpha), 
                          shape=shape, color=col, size = pointsize)
    else if("arrow" %in% geom) 
      p <- p + geom_segment(data = data,
                            aes_string(x = 0, y = 0, xend = x, yend = y, alpha=alpha),
                            arrow = grid::arrow(length = grid::unit(0.2, 'cm')), color=col)
    
    if(lab & "text" %in% geom) 
      p <- p + geom_text(data = label_coord, 
                         aes_string(x, y, alpha=alpha, label="name"),
                         size = labelsize, color=col, vjust = -0.7)
    
    if(!is.null(alpha.limits)) p <- p + scale_alpha(limits = alpha.limits)
  }
  
  else{
    if("point" %in% geom) 
      p <- p + geom_point(data = data, aes_string(x, y),
                          shape=shape, color=col, size = pointsize)
    else if("arrow" %in% geom) 
      p <- p + geom_segment(data = data,
                            aes_string(x = 0, y = 0, xend = x, yend = y),
                            arrow = grid::arrow(length = grid::unit(0.2, 'cm')), color=col)
    if(lab & "text" %in% geom) 
      p <- p + geom_text(data = label_coord, aes_string(x,y), 
                         color = col, label = data$name, size = labelsize, vjust = -0.7)
  }
  
  return(p)
}

# Make arrows from the origine to the point (x, y)
#+++++++++++++++++++++++++++++++++++
## p: a ggplot. If not null, points are added on the existing plot
## data: is a data frame of form: name|x|y|coord|cos2|contrib
## col: point color. The default value is "black".
#  Possible values include also : "cos2", "contrib", "coord", "x" or "y".
#  In this case, the colors for variables are automatically controlled by their qualities ("cos2"),
#  contributions ("contrib"), coordinates (x^2+y2, "coord"), x values("x") or y values("y")
## x, y: contains the name of x and y variables
## alpha: controls the transparency of colors.
# The value can variate from 0 (total transparency) to 1 (no transparency).
# Default value is 1. Possible values include also : "cos2", "contrib", "coord", "x" or "y".
## alpha.limits: the range of alpha values. Used only when alpha is 
# a continuous variable(i.e, cos2, contrib, ...)
## lab: if TRUE points are labelled
## geom: a character specifying the geometry to be used for the graph.
#  Allowed values are "point" (to show only points),
#  "text" to show only labels or c("point", "text") to show both types.
.ggarrow <- function(p = NULL, data, x = 'x', y = 'y', col="black",  alpha = 1,
                       alpha.limits = NULL,
                       shape = "19", 
                       geom=c("arrow", "text"), lab = TRUE, labelsize = 4,
                       jitter = list(what = "label", width = NULL, height = NULL),... )
{
  
  data <- as.data.frame(data)
  label_coord <- data
  
  # jittering
  if(jitter$what %in% c("both", "b")){
    label_coord <- data <- .jitter(data, jitter)
  }
  else if(jitter$what %in% c("point", "p")){
    data <- .jitter(data, jitter)
  }
  else if(jitter$what %in% c("label", "l")){
    label_coord <- .jitter(label_coord, jitter)
  }
  
  if(is.null(p)) p <- ggplot() 
  # The color and the transparency of variables are automatically controlled by
  # their cos2, contrib,  "x" or "y" coordinates
  if(col %in% c("cos2","contrib", "coord", "x", "y") &
       alpha %in% c("cos2","contrib", "coord", "x", "y"))
  {
    if("point" %in% geom) 
      p <- p + geom_point(data = data, 
                          aes_string(x,y, color=col, alpha=alpha), shape=shape)
    
    if(lab & "text"%in% geom) 
      p <- p + geom_text(data = label_coord, 
                         aes_string(x,y, label = 'name', 
                                    color=col, alpha=alpha), size = labelsize)
    if(!is.null(alpha.limits)) p <- p + scale_alpha(limits = alpha.limits)
  }
  # Only the color is controlled automatically
  else if(col %in% c("cos2","contrib", "coord", "x", "y")){
    
    if("point" %in% geom) 
      p <- p + geom_point(data = data, aes_string(x,y, color=col),
                          shape=shape, alpha=alpha)
    
    if(lab & "text" %in% geom) 
      p <- p + geom_text(data = label_coord,
                         aes_string(x, y, color=col),
                         label = data$name,  size = labelsize, alpha=alpha, vjust = -0.7)
  }
  
  # Only the transparency is controlled automatically
  else if(alpha %in% c("cos2","contrib", "coord", "x", "y")){
    
    if("point" %in% geom) 
      p <- p + geom_point(data = data, 
                          aes_string(x, y, alpha=alpha), shape=shape, color=col)
    
    if(lab & "text" %in% geom) 
      p <- p + geom_text(data = label_coord, 
                         aes_string(x, y, alpha=alpha, label="name"),
                         size = labelsize, color=col, vjust = -0.7)
    
    if(!is.null(alpha.limits)) p <- p + scale_alpha(limits = alpha.limits)
  }
  
  else{
    if("point" %in% geom) 
      p <- p + geom_point(data = data, aes(x, y), shape=shape, color=col)
    if(lab & "text" %in% geom) 
      p <- p + geom_text(data = label_coord, aes_string(x,y), 
                         color = col, label = data$name, size = labelsize, vjust = -0.7)
  }
  
  return(p)
}



# Points jitter, to reduce overploting
# data : a result from facto_summarize containing the x and y coordinates of points
# jitter: a vector of length 3 containing the width and the height of jitter and the 
# element to be jittered ("label", "point", "both")
.jitter <- function(data, jitter = list(what = "label", width = NULL, height = NULL)){
  
  if(!is.null(data)){
    if(!is.null(jitter$width)){
      width <- abs(jitter$width)
     set.seed(1234)
      xjit <- runif(nrow(data), min = -width, max = width)
      data$x <- data$x + ((xjit + rev(xjit))/2)
    }
    
    if(!is.null(jitter$height)){
      height <- abs(jitter$height)
       set.seed(12345)
      yjit <- runif(nrow(data), min = -height, max = height)
      data$y <- data$y + ((yjit + rev(yjit))/2)
    }
    
  }
  
  return(data)
}



# Finih a plot
#++++++++++++++++++++++++++++++
# p a ggplot2
# X an object of class PCA, MCA or CA
# axes the plotted axes
.fviz_finish <- function(p, X, axes = 1:2){
  
  cc <- .get_facto_class(X)
  title <- paste0(cc, " factor map")
  
  eig <- get_eigenvalue(X)[,2]
  xlab = paste0("Dim", axes[1], " (", round(eig[axes[1]],1), "%)") 
  ylab = paste0("Dim", axes[2], " (", round(eig[axes[2]], 1),"%)")
  
  p <- p +
    geom_hline(yintercept = 0, color = "black", linetype="dashed") +
    geom_vline(xintercept = 0, color = "black", linetype="dashed") +
    labs(title = title, x = xlab, y = ylab)
  
  return(p)
}

# Check the element to be labelled
#+++++++++++++++++
## label: a character vector specifying the elements to be labelled
# possible values are "all", "none" or the combination of 
# c("row", "row.sup", "col", "col.sup", "ind", "ind.sup", "quali", "var", "quanti.sup")
## Returns a list 
.label <- function(label){
  lab  <- list()
  # var - PCA, MCA
  lab$var <- lab$quanti <- lab$quali.sup <- FALSE
  if(label[1]=="all" | "var" %in% label) lab$var =TRUE
  if(label[1]=="all" | "quanti.sup" %in% label) lab$quanti =TRUE
  if(label[1]=="all" | "quali.sup" %in% label) lab$quali.sup =TRUE
  # ind - PCA, MCA
  lab$ind <- lab$ind.sup <- lab$quali <- FALSE
  if(label[1]=="all" | "ind" %in% label) lab$ind =TRUE
  if(label[1]=="all" | "ind.sup" %in% label) lab$ind.sup =TRUE
  if(label[1]=="all" | "quali" %in% label) lab$quali =TRUE
  # row - ca
  lab$row <- lab$row.sup <- FALSE
  if(label[1]=="all" | "row" %in% label) lab$row =TRUE
  if(label[1]=="all" | "row.sup" %in% label) lab$row.sup =TRUE
  # col - ca
  lab$col <- lab$col.sup <- FALSE
  if(label[1]=="all" | "col" %in% label) lab$col =TRUE
  if(label[1]=="all" | "col.sup" %in% label) lab$col.sup =TRUE
  
  lab
}


# Check the element to be hidden
#+++++++++++++++++
## invisible: a character vector specifying the elements to be hidden
# possible values are "all", "none" or the combination of 
# c("row", "row.sup", "col", "col.sup", "ind", "ind.sup", "quali", "quali.sup",  "var", "quanti.sup")
## Returns a list 
.hide <- function(invisible){
  hide  <- list()
  # var - PCA, MCA
  hide$var <- hide$quanti <- hide$quali.sup <- FALSE
  if("var" %in% invisible) hide$var =TRUE
  if("quanti.sup" %in% invisible) hide$quanti =TRUE
  if("quali.sup" %in% invisible) hide$quali.sup = TRUE
  # ind - PCA, MCA
  hide$ind <- hide$ind.sup <- hide$quali <- FALSE
  if("ind" %in% invisible) hide$ind =TRUE
  if("ind.sup" %in% invisible) hide$ind.sup =TRUE
  if("quali" %in% invisible) hide$quali =TRUE
  # row - ca
  hide$row <- hide$row.sup <- FALSE
  if("row" %in% invisible) hide$row =TRUE
  if("row.sup" %in% invisible) hide$row.sup =TRUE
  # col - ca
  hide$col <- hide$col.sup <- FALSE
  if("col" %in% invisible) hide$col =TRUE
  if("col.sup" %in% invisible) hide$col.sup =TRUE
  
  hide
}





