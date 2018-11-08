library( lattice )

### Preparing data

NHANESnames <- c( "BIOPRO_G.XPT", "CBC_G.XPT", "DEMO_G.XPT", "GHB_G.XPT", "HDL_G.XPT" )
BaseURL <- "https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/"

nhanes <- Reduce( function( ... ) merge( ... ),
                  lapply( NHANESnames, function( n ) Hmisc::sasxport.get( paste0( BaseURL, n ), lowernames = FALSE ) ) )

nhanes <- nhanes[ !is.na( nhanes$RIDAGEYR )&nhanes$RIDAGEYR>18, ]
nhanes <- nhanes[ !is.na( nhanes$RIAGENDR )&nhanes$RIAGENDR==1, ]

nhanes <- data.frame( ANC = nhanes$LBDNENO, ABC = nhanes$LBDBANO,
                      ALC = nhanes$LBDLYMNO, AMC = nhanes$LBDMONO,
                      AEC = nhanes$LBDEONO, RBC = nhanes$LBXRBCSI,
                      HGB = nhanes$LBXHGB * 10, HCT = nhanes$LBXHCT / 100,
                      MCV = nhanes$LBXMCVSI, MCH = nhanes$LBXMCHSI,
                      MCHC = nhanes$LBXMC * 10, RDW = nhanes$LBXRDW,
                      PLT = nhanes$LBXPLTSI, MPV = nhanes$LBXMPSI,
                      SNA = nhanes$LBXSNASI, GHB = nhanes$LBXGH / 100,
                      SK = nhanes$LBXSKSI, SCL = nhanes$LBXSCLSI,
                      SCA = nhanes$LBDSCASI, SP = nhanes$LBDSPHSI,
                      CPK = nhanes$LBXSCK, STB = nhanes$LBDSTBSI,
                      BIC = nhanes$LBXSC3SI, GLU = nhanes$LBDSGLSI,
                      IRN = nhanes$LBDSIRSI, LDH = nhanes$LBXSLDSI,
                      STP = nhanes$LBDSTPSI, SUA = nhanes$LBDSUASI,
                      SAL = nhanes$LBDSALSI, TRI = nhanes$LBDSTRSI,
                      SGL = nhanes$LBDSGBSI, BUN = nhanes$LBDSBUSI,
                      SCR = nhanes$LBDSCRSI, STC = nhanes$LBDSCHSI,
                      HDL = nhanes$LBDHDDSI, AST = nhanes$LBXSASSI,
                      ALT = nhanes$LBXSATSI, GGT = nhanes$LBXSGTSI,
                      ALP = nhanes$LBXSAPSI, WEIGHT = nhanes$WTMEC2YR )

nhanes <- na.omit( nhanes )

nhanesWeight <- nhanes$WEIGHT
nhanes <- nhanes[ , -40 ]

###  Preliminary calculations

p <- ncol( nhanes )

Binned <- nhanes
for( i in 1:p ) {
  Binned[ , i ] <- cut( Binned[ , i ], breaks = 4 )
}

TotalCorrelation <- function( RawData, Weights, IndEntropies = NULL ) {
  if( is.null( IndEntropies ) ) {
    IndEntropies <- apply( RawData, 2, function( x ) {
      entropy::entropy.shrink( survey::svytable( ~ . , survey::svydesign( ~ 0, data = as.data.frame( x ), weights = Weights ) ),
                               verbose = FALSE ) } )
  }
  sum( IndEntropies ) - entropy::entropy.shrink( survey::svytable( ~ . , survey::svydesign( ~ 0, data = RawData,
                                                                                            weights = Weights ) ) )
}

PearsonCorMat <- cov.wt( nhanes, nhanesWeight, cor = TRUE )$cor
TotCorMat <- matrix( apply( expand.grid( x = 1:p, y = 1:p ), 1,
                            function( x ){ return( TotalCorrelation( Binned[ , c( x[ 1 ], x[ 2 ] ) ], nhanesWeight ) ) } ),
                     nc = p )
colnames( TotCorMat ) <- rownames( TotCorMat ) <- colnames( PearsonCorMat )

### Visualizations

range( TotCorMat )

levelplot( PearsonCorMat, col.regions = heat.colors( 100 ), xlab = "", ylab = "", scales = list( x = list( rot = 90 ) ) )
levelplot( PearsonCorMat, col.regions = gray( 0:100/100 ), xlab = "", ylab = "", scales = list( x = list( rot = 90 ) ) )
levelplot( abs( PearsonCorMat ), col.regions = heat.colors( 100 ), xlab = "", ylab = "", scales = list( x = list( rot = 90 ) ) )
cairo_pdf( "PearsonCorr.pdf" )
levelplot( abs( PearsonCorMat ), col.regions = gray( 0:100/100 ), colorkey = TRUE,
           at = exp( seq( log( 1e-05 ), log( 1 ), length.out = 100 ) ),
           xlab = "", ylab = "", scales = list( x = list( rot = 90 ) ) )
dev.off()

levelplot( TotCorMat, col.regions = heat.colors( 100 ), xlab = "", ylab = "", scales = list( x = list( rot = 90 ) ) )
levelplot( log( TotCorMat ), col.regions = heat.colors( 100 ), xlab = "", ylab = "", scales = list( x = list( rot = 90 ) ) )
levelplot( TotCorMat, col.regions = heat.colors( 100 ), colorkey = TRUE, scales = list( x = list( rot = 90 ) ),
           at = exp( seq( log( 1e-07 ), log( 1 ), length.out = 100 ) ), xlab = "", ylab = "" )
cairo_pdf( "TotCorrMat.pdf" )
levelplot( TotCorMat, col.regions = gray( 0:100/100 ), colorkey = TRUE,
           at = exp( seq( log( 1e-07 ), log( 1 ), length.out = 100 ) ),
           scales = list( x = list( rot = 90 ) ), xlab = "", ylab = "" )
dev.off()

### First approach: plain old hierachical clustering

cairo_pdf( "PearsonDendr.pdf" )
plot( as.dendrogram( hclust( as.dist( 1 - log( abs( PearsonCorMat ) ) ), method = "ward.D2" ), hang = 0.1 ),
      main = "", sub = "", xlab = "", ylim = c( 0, 16 ) )
dev.off()

cairo_pdf( "TotCorrDendr.pdf" )
plot( as.dendrogram( hclust( as.dist( 1 - log( TotCorMat ) ), method = "ward.D2" ), hang = 0.1 ),
      main = "", sub = "", xlab = "", ylim = c( 0, 16 ) )
dev.off()

### Second approach: greedy growth

DecodeLUT <- function( LUT, Code, p ) {
  if( Code <= p ) {
    return( Code )
  } else {
    return( c( DecodeLUT( LUT, LUT[ LUT[ , 1 ] == Code, 2 ], p ), DecodeLUT( LUT, LUT[ LUT[ , 1 ] == Code, 3 ], p ) ) )
  }
}

ToMerge <- matrix( rep( 1:p, p ), nc = p, byrow = TRUE )
LookupNew <- matrix( c( seq( p + 1, 2 * p ), rep( 0, 3 * p ) ), nc = 4 )

MatGrid <- t( combn( 1:p, 2 ) )
TotCorMat <- apply( MatGrid , 1, function( x ){ return( TotalCorrelation( Binned[ , c( x[ 1 ], x[ 2 ] ) ], nhanesWeight ) ) } )

hist( log10( TotCorMat ), xlab = "log10(Total correlation)", main = "" )

MaxPair <- as.numeric( MatGrid[ which.max( TotCorMat ), ] )

ToMerge[ 2, MaxPair ] <- p + 1
LookupNew[ 1, 2:3 ] <- MaxPair
LookupNew[ 1, 4 ] <- max( TotCorMat )

for ( i in 2:( p - 1 ) ) {
  print(i)
  MatGrid <- expand.grid( x = unique( ToMerge[ i, ] ), y = unique( ToMerge[ i, ] ) )
  MatGrid <- MatGrid[ MatGrid$x > MatGrid$y, ]
  TotCorMat <- as.numeric( apply( MatGrid, 1, function( x ){ return( TotalCorrelation(
    Binned[ , c( DecodeLUT( LookupNew, x[ 1 ], p ), DecodeLUT( LookupNew, x[ 2 ], p ) ) ], nhanesWeight ) ) } ) )
  MaxPair <- as.numeric( MatGrid[ which.max( TotCorMat ), ] )
  ToMerge[ i + 1, ] <- ToMerge[ i, ]
  ToMerge[ i + 1, c( DecodeLUT( LookupNew, MaxPair[ 1 ], p ), DecodeLUT( LookupNew, MaxPair[ 2 ], p ) ) ] <- p + i
  LookupNew[ i, 2:3 ] <- MaxPair
  LookupNew[ i, 4 ] <- max( TotCorMat )
}

sapply( apply( LookupNew[ 1:11, 2:3 ], 1, function( x ){ return( c( DecodeLUT( LookupNew, x[ 1 ], p ),
                                                                    DecodeLUT( LookupNew, x[ 2 ], p ) ) ) } ),
        function( x ){ return( names( nhanes[ x ] ) ) } )

# Third approach: APRIORI-style search tree, with pruning

EntropyOnedim <- apply( Binned, 2, function( x ) {
  entropy::entropy.shrink( survey::svytable( ~ . , survey::svydesign( ~ 0, data = as.data.frame( x ),
                                                                      weights = nhanesWeight ) ), verbose = FALSE ) } )
nms <- names( nhanes )
Limit <- 0.1

MatGrid <- t( combn( nms, 2 ) )
TotCorMat <- data.frame( MatGrid, TotCor = apply( MatGrid, 1, function( x ){
  TotalCorrelation( Binned[ , c( x[ 1 ], x[ 2 ] ) ], nhanesWeight, EntropyOnedim[ c( x[ 1 ], x[ 2 ] ) ] ) } ),
  stringsAsFactors = FALSE )
cairo_pdf( "HistPairTotCorr.pdf" )
hist( log10( TotCorMat$TotCor ), xlab = "log10(Total correlation)", main = "" )
dev.off()
TotCorMat <- TotCorMat[ ( TotCorMat$TotCor > Limit ), ]
TotCorMat

nms <- setdiff( nms, c( TotCorMat$X1, TotCorMat$X2 ) )
MatGrid <- t( combn( nms, 3 ) )
TotCorMat <- data.frame( MatGrid, TotCor = apply( MatGrid, 1, function( x ){
  TotalCorrelation( Binned[ , c( x[ 1 ], x[ 2 ], x[ 3 ] ) ], nhanesWeight, EntropyOnedim[ c( x[ 1 ], x[ 2 ], x[ 3 ] ) ] ) } ),
  stringsAsFactors = FALSE )
hist( log10( TotCorMat$TotCor ), xlab = "log10(Total correlation)", main = "" )
TotCorMat <- TotCorMat[ ( TotCorMat$TotCor > Limit ), ]
TotCorMat

nms <- setdiff( nms, c( TotCorMat$X1, TotCorMat$X2, TotCorMat$X3 ) )
MatGrid <- t( combn( nms, 4 ) )
TotCorMat <- data.frame( MatGrid, TotCor = apply( MatGrid, 1, function( x ){
  TotalCorrelation( Binned[ , c( x[ 1 ], x[ 2 ], x[ 3 ], x[ 4 ] ) ], nhanesWeight,
                    EntropyOnedim[ c( x[ 1 ], x[ 2 ], x[ 3 ], x[ 4 ] ) ] ) } ), stringsAsFactors = FALSE )
hist( log10( TotCorMat$TotCor ), xlab = "log10(Total correlation)", main = "" )
TotCorMat <- TotCorMat[ ( TotCorMat$TotCor > Limit ), ]
TotCorMat

nms <- setdiff( nms, c( TotCorMat$X1, TotCorMat$X2, TotCorMat$X3, TotCorMat$X4 ) )
MatGrid <- t( combn( nms, 5 ) )
TotCorMat <- data.frame( MatGrid, TotCor = apply( MatGrid, 1, function( x ){
  TotalCorrelation( Binned[ , c( x[ 1 ], x[ 2 ], x[ 3 ], x[ 4 ], x[ 5 ] ) ], nhanesWeight,
                    EntropyOnedim[ c( x[ 1 ], x[ 2 ], x[ 3 ], x[ 4 ], x[ 5 ] ) ] ) } ), stringsAsFactors = FALSE )
hist( log10( TotCorMat$TotCor ), xlab = "log10(Total correlation)", main = "" )
TotCorMat <- TotCorMat[ ( TotCorMat$TotCor > Limit ), ]
TotCorMat
