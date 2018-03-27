#==================================================
#==================================================
#==================================================
# Libraries
require( bio3d )
library( foreach )
library( doParallel )



#==================================================
#==================================================
#==================================================
# Functions

#==================================================
# Check if is letter
is.letter <- function( x )
{
  grepl( "[[:alpha:]]",
         x )
}

# #==================================================
# # Main CRISPro function
# CRISPro <- function( overlap.cutoff = 1.5,
#                      distance_angstroem = 10,
#                      cutoff = 0.0005,
#                      file.rama = NULL,
#                      file.exe.mutate = NULL,
#                      output_prefix = "OUTPUT",
#                      option = NULL,
#                      file.pdb.prefusion = NULL,
#                      prefusion.chain = NULL,
#                      file.pdb.postfusion = NULL,
#                      postfusion.chain = NULL )
# {
#   
# }

#==================================================
#==================================================
#==================================================
# Main

#==================================================
# start time measurement
start.main <- proc.time()
#==================================================


# #==================================================
# # Set working directory
# setwd( "~/Documents/work/vrc/sbis_projects/CRISPro/data/RSVF/option-1" )

#==================================================
# Source function

#==================================================
# Command line arguments
overlap.cutoff <- as.numeric( commandArgs()[ 3 ] )
distance_angstroem <- as.numeric( commandArgs()[ 4 ] )
cutoff <- as.numeric( commandArgs()[ 5 ] )
n.cores <- as.numeric( commandArgs()[ 6 ] )
registerDoParallel( n.cores )
file.rama <- commandArgs()[ 7 ]
file.exe.mutate <- commandArgs()[ 8 ]
output_prefix <- commandArgs()[ 9 ]
option <- as.numeric( commandArgs()[ 10 ] )
if( option == 1 )
{
  file.pdb.prefusion <- commandArgs()[ 11 ]
  prefusion.chain <- commandArgs()[ 12 ]
  file.pdb.postfusion <- commandArgs()[ 13 ]
  postfusion.chain <- commandArgs()[ 14 ]
} else if( option == 2 )
{
  file.pdb.prefusion <- commandArgs()[ 11 ]
  prefusion.chain <- commandArgs()[ 12 ]
} else if( option == 3 )
{
  file.pdb.postfusion <- commandArgs()[ 11 ]
  postfusion.chain <- commandArgs()[ 12 ]
} else
{
  stop( paste( "ERROR: Option", option, "is not an input. Please choose 1, 2, or 3!" ) )
}

# overlap.cutoff <- 1.5
# distance_angstroem <- 10
# cutoff <- 0.0005
# n.cores <- 8
# registerDoParallel( n.cores )
# file.rama <- "~/Documents/work/vrc/sbis_projects/CRISPro/scripts/rama8000-transpro.data"
# file.exe.mutate <- "~/Documents/work/vrc/sbis_projects/CRISPro/scripts/mutate.py"
# output_prefix <- "RSVF_option-1"
# option <- 1
# if( option == 1 )
# {
#   file.pdb.prefusion <- "4jhw_trimer.pdb"
#   prefusion.chain <- "A"
#   file.pdb.postfusion <- "3rrr_trimer.pdb"
#   postfusion.chain <- "A"
# } else if( option == 2 )
# {
#   file.pdb.prefusion <- commandArgs()[ 11 ]
#   prefusion.chain <- commandArgs()[ 12 ]
# } else if( option == 3 )
# {
#   file.pdb.postfusion <- commandArgs()[ 11 ]
#   postfusion.chain <- commandArgs()[ 12 ]
# } else
# {
#   stop( paste( "ERROR: Option", option, "is not an input. Please choose 1, 2, or 3!" ) )
# }


#==================================================
# Run CRISPro

#==================================================
# Set variables
file.trans.proline <- file.rama
vec.potential.x <- seq( -179,
                        179,
                        2 )


#==================================================
# Load trans proline phi psi angle lattice
df.trans.proline <- data.frame( read.table( file.trans.proline ) )
colnames( df.trans.proline ) <- c( "x",
                                   "y",
                                   "density" )
x <- df.trans.proline[ df.trans.proline$density > cutoff, ]$x
y <- df.trans.proline[ df.trans.proline$density > cutoff, ]$y
df.x.y <- data.frame( cbind( x,
                             y ) )


#==================================================
#==================================================
#==================================================
# Prefusion conformation experiments
if( option == 1 )
{
  #==================================================
  #==================================================
  #==================================================
  # Prefusion conformation
  
  #==================================================
  # Load original PDB and output altered one removing all non-ATOM entries
  
  # Load input PDB
  pdb.prefusion.raw <- read.pdb( file.pdb.prefusion )
  pdb.prefusion.raw.atom <- pdb.prefusion.raw$atom
  
  # Save new altered PDB
  pdb.prefusion.raw.atom.altered <- pdb.prefusion.raw.atom[ pdb.prefusion.raw.atom$type == "ATOM", ]
  file.pdb.prefusion.altered <- paste( output_prefix,
                                       "_prefusion.pdb",
                                       sep = "" )
  write.pdb( file = file.pdb.prefusion.altered,
             xyz = as.vector( t( cbind( pdb.prefusion.raw.atom.altered$x,
                                        pdb.prefusion.raw.atom.altered$y,
                                        pdb.prefusion.raw.atom.altered$z ) ) ),
             type = pdb.prefusion.raw.atom.altered$type,
             eleno = pdb.prefusion.raw.atom.altered$eleno,
             elety = pdb.prefusion.raw.atom.altered$elety,
             alt = pdb.prefusion.raw.atom.altered$alt,
             resid = pdb.prefusion.raw.atom.altered$resid,
             chain = pdb.prefusion.raw.atom.altered$chain,
             resno = pdb.prefusion.raw.atom.altered$resno,
             insert = pdb.prefusion.raw.atom.altered$insert,
             o = pdb.prefusion.raw.atom.altered$o,
             b = pdb.prefusion.raw.atom.altered$b )
  
  pdb.prefusion <- read.pdb( file.pdb.prefusion.altered )
  pdb.prefusion.atom <- pdb.prefusion$atom
  
  #==================================================
  # Generate data frame that saves chain and resno
  df.chain.resno.resid.prefusion <- NULL
  for( j in 1:nrow( pdb.prefusion.atom ) )
  {
    var.chain <- pdb.prefusion.atom[ j, ]$chain
    if( is.na( pdb.prefusion.atom[ j, ]$insert ) )
    {
      var.resno <- as.character( pdb.prefusion.atom[ j, ]$resno )
    } else
    {
      var.resno <- paste( pdb.prefusion.atom[ j, ]$resno,
                          pdb.prefusion.atom[ j, ]$insert,
                          sep = "" )
    }
    var.resid <- pdb.prefusion.atom[ j, ]$resid
    
    if( length( df.chain.resno.resid.prefusion[ df.chain.resno.resid.prefusion[ , 1 ] == var.chain &
                                                df.chain.resno.resid.prefusion[ , 2 ] == var.resno, ] ) == 0 )
    {
      df.chain.resno.resid.prefusion <- rbind( df.chain.resno.resid.prefusion,
                                               c( var.chain,
                                                  var.resno,
                                                  var.resid ) )
    }
  }
  df.chain.resno.resid.prefusion <- data.frame( df.chain.resno.resid.prefusion )
  colnames( df.chain.resno.resid.prefusion ) <- c( "chain_prefusion",
                                                   "resno_prefusion",
                                                   "resid_prefusion" )
  df.chain.resno.resid.prefusion$chain_prefusion <- as.vector( df.chain.resno.resid.prefusion$chain_prefusion )
  df.chain.resno.resid.prefusion$resno_prefusion <- as.vector( df.chain.resno.resid.prefusion$resno_prefusion )
  df.chain.resno.resid.prefusion$resid_prefusion <- as.vector( aa321( df.chain.resno.resid.prefusion$resid_prefusion ) )
  
  #==================================================
  # Calculate torsion angles for PDB using DSSP
  
  # Determine torsion angles (and secondary structure) using DSSP
  list.ss.dssp.prefusion <- dssp( pdb.prefusion,
                                  resno = TRUE ) 
  vec.phi.prefusion <- list.ss.dssp.prefusion$phi
  vec.phi.prefusion[ which( vec.phi.prefusion == 360 ) ] <- NA
  
  vec.psi.prefusion <- list.ss.dssp.prefusion$psi
  vec.psi.prefusion[ which( vec.psi.prefusion == 360 ) ] <- NA
  
  # Select torsion angles of interest
  df.prefusion.phi.psi <- data.frame( cbind( vec.phi.prefusion,
                                             vec.psi.prefusion ) )
  
  #==================================================
  # Combine PDB information with torsion angles
  df.chain.resno.resid.prefusion.phi.psi <- data.frame( cbind( df.chain.resno.resid.prefusion,
                                                               df.prefusion.phi.psi ) )
  colnames( df.chain.resno.resid.prefusion.phi.psi ) <- c( colnames( df.chain.resno.resid.prefusion ),
                                                           "phi_prefusion",
                                                           "psi_prefusion" )
  
  #==================================================
  # Identify if angles are in trans proline conformation
  vec.phi.psi.prefusion.trans.proline <- NULL
  for( i in 1:nrow( df.chain.resno.resid.prefusion.phi.psi) )
  {
    phi <- df.chain.resno.resid.prefusion.phi.psi[ i, ]$phi_prefusion
    psi <- df.chain.resno.resid.prefusion.phi.psi[ i, ]$psi_prefusion
    
    if( is.na( phi ) | is.na( psi ) )
    {
      vec.phi.psi.prefusion.trans.proline <- c( vec.phi.psi.prefusion.trans.proline,
                                                NA )
    } else
    {
      phi.tmp <- abs( phi - vec.potential.x )
      phi.tmp.choice <- which( phi.tmp == min( phi.tmp ) )
      if( length( phi.tmp.choice ) > 1 )
      {
        # phi.tmp.choice <- sample( phi.tmp.choice,
        #                           1 )
        phi.tmp.choice <- min( phi.tmp.choice )
      }
      phi.round <- vec.potential.x[ phi.tmp.choice ]
      phi.round
      
      psi.tmp <- abs( psi - vec.potential.x )
      psi.tmp.choice <- which( psi.tmp == min( psi.tmp ) )
      if( length( psi.tmp.choice ) > 1 )
      {
        # psi.tmp.choice <- sample( psi.tmp.choice,
        #                           1 )
        psi.tmp.choice <- min( psi.tmp.choice )
      }
      psi.round <- vec.potential.x[ psi.tmp.choice ]
      psi.round
      
      if( nrow( df.x.y[ df.x.y$x == phi.round &
                        df.x.y$y == psi.round, ] ) > 0 )
      {
        vec.phi.psi.prefusion.trans.proline <- c( vec.phi.psi.prefusion.trans.proline,
                                                  1 )
      } else
      {
        vec.phi.psi.prefusion.trans.proline <- c( vec.phi.psi.prefusion.trans.proline,
                                                  0 )
      }
    }
  }
  
  # Add column to data frame
  df.chain.resno.resid.prefusion.phi.psi.angle <- cbind( df.chain.resno.resid.prefusion.phi.psi,
                                                         vec.phi.psi.prefusion.trans.proline )
  colnames( df.chain.resno.resid.prefusion.phi.psi.angle ) <- c( colnames( df.chain.resno.resid.prefusion.phi.psi ),
                                                                 "angle_trans_proline_prefusion" )
  
  
  
  #==================================================
  # Secondary structure and PHI and PSI angles
  
  # Add extra column to add helix conformations
  df.chain.resno.resid.prefusion.phi.psi.angle.helix <- cbind( df.chain.resno.resid.prefusion.phi.psi.angle,
                                                               rep( 0,
                                                                    nrow( df.chain.resno.resid.prefusion.phi.psi.angle ) ),
                                                               rep( 0,
                                                                    nrow( df.chain.resno.resid.prefusion.phi.psi.angle ) ),
                                                               rep( 0,
                                                                    nrow( df.chain.resno.resid.prefusion.phi.psi.angle ) ) )
  colnames( df.chain.resno.resid.prefusion.phi.psi.angle.helix ) <- c( colnames( df.chain.resno.resid.prefusion.phi.psi.angle ),
                                                                       "helix_prefusion",
                                                                       "helix_termini_prefusion",
                                                                       "DSSP_prefusion" )
  
  
  # Select subsets where chain of interest if given
  vec.ind.subset.chain <- as.numeric( which( list.ss.dssp.prefusion$helix$chain == prefusion.chain ) )
  
  vec.helix.start <- list.ss.dssp.prefusion$helix$start[ vec.ind.subset.chain ]
  vec.helix.end <- list.ss.dssp.prefusion$helix$end[ vec.ind.subset.chain ]
  vec.index.helix <- NULL
  if( length( vec.helix.start ) > 0 )
  {
    for( i in 1:length( vec.helix.start ) )
    {
      vec.index.helix <- c( vec.index.helix,
                            vec.helix.start[ i ]:vec.helix.end[ i ] )
    }
    
    vec.resno.helix <- df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% as.character( vec.index.helix ), ]$resno_prefusion
    
    #==================================================
    # Add to data frame
    df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% vec.resno.helix, ]$helix_prefusion <- 1
    # df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% c( vec.helix.start, vec.helix.end ), ]$helix_termini_prefusion <- 1
    df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% as.character( sort( as.vector( c( vec.helix.start, vec.helix.end ) ) ) ), ]$helix_termini_prefusion <- 1
  }
  
  # Add DSSP secondary structure to data frame
  vec.dssp <- list.ss.dssp.prefusion$sse
  vec.dssp[ vec.dssp == " " ] = "L"
  df.chain.resno.resid.prefusion.phi.psi.angle.helix$DSSP_prefusion <- vec.dssp
  
  #==================================================
  # Run pymol to mutate all positions
  print( "Run pymol to mutate all positions in the prefusion conformation" )
  df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain <- df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$chain_prefusion == prefusion.chain, ]
  
  # vec.prefusion.no.clashes <- NULL
  # for( i in 1:nrow( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ) )
  # {
  #   start.time <- proc.time()
  #   print( paste( i,
  #                 nrow( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ),
  #                 sep = "/" ) )
  vec.prefusion.no.clashes <- foreach( i = 1:nrow( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ), .combine = c ) %dopar%
  {
    chain <- df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain[ i, 1 ]
    resno <- df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain[ i, 2 ]
    selection <- paste( chain,
                        "/",
                        resno,
                        "/",
                        sep = "" )
    
    file.pdb.prefusion.mutation <- paste( unlist( strsplit( file.pdb.prefusion,
                                                            "\\.pdb" ) ),
                                          "_",
                                          chain,
                                          resno,
                                          ".pdb",
                                          sep = "" )
    
    #==================================================
    # Generate mutated PDB file
    command <- paste( "pymol",
                      "-qc",
                      # "~/Documents/work/vrc/sbis_projects/proline_predictor/scripts/mutate.py",
                      # "mutate.py",
                      file.exe.mutate,
                      "--",
                      file.pdb.prefusion,
                      file.pdb.prefusion.mutation,
                      selection,
                      "PRO",
                      distance_angstroem )
    # Run python (PyMol) command
    system( command,
            intern = TRUE )
    
    #==================================================
    # Analyze mutated protein structure
    
    # Load mutated structure
    pdb.mutant <- read.pdb( file.pdb.prefusion.mutation )
    pdb.mutant.atom <- pdb.mutant$atom
    
    # Select mutated residue
    # Do we have an insertion?
    if( is.letter( resno ) )
    {
      resno.new <- NULL
      insert.new <- NULL
      for( char in unlist( strsplit( resno,
                                     "" ) ) )
      {
        if( is.letter( char ) )
        {
          insert.new <- c( insert.new,
                           char )
        } else
        {
          resno.new <- c( resno.new,
                          char )
        }
      }
      
      resno.new <- as.numeric( paste( resno.new,
                                      collapse = "" ) )
      insert.new <- paste( insert.new,
                           collapse = "" )
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno.new &
                                                     pdb.mutant.atom$insert == insert.new &
                                                     !( is.na( pdb.mutant.atom$insert ) ), ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno.new &
                                                            pdb.mutant.atom$insert == insert.new &
                                                            !( is.na( pdb.mutant.atom$insert ) ) ), ]
      
    } else
    {
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno, ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno ), ]
    }
    # Select xyz coordinates for non-backbone atoms of mutated residue
    pdb.mutant.atom.mutation.no.backbone <- pdb.mutant.atom.mutation[ which( !( pdb.mutant.atom.mutation$elety %in% c( "N",
                                                                                                                       "CA",
                                                                                                                       "C",
                                                                                                                       "O" ) ) ), ]
    pdb.mutant.atom.mutation.no.backbone.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.mutation.no.backbone[ , 9:11 ] ) ) )
    
    # Select xyz coordinates for all other residues
    pdb.mutant.atom.non.mutation.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.non.mutation[ , 9:11 ] ) ) )
    
    # Calculate distance matrix
    mat.dist <- dist.xyz( pdb.mutant.atom.mutation.no.backbone.xyz,
                          pdb.mutant.atom.non.mutation.xyz )
    
    # Generate vectors representing atom types for mutated residue and rest
    vec.atom.type.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.mutation.no.backbone$elety ) )
    {
      elety <- pdb.mutant.atom.mutation.no.backbone$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.mutation <- c( vec.atom.type.mutation,
                                   elety.split[ which( unlist( lapply( elety.split,
                                                                       is.letter ) ) )[ 1 ] ] )
    }
    
    vec.atom.type.non.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.non.mutation$elety ) )
    {
      elety <- pdb.mutant.atom.non.mutation$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.non.mutation <- c( vec.atom.type.non.mutation,
                                       elety.split[ which( unlist( lapply( elety.split,
                                                                           is.letter ) ) )[ 1 ] ] )
    }
    
    # Calculate MolProbity clash constraint
    df.atom.type.mutation.nonmutation <- data.frame( expand.grid( unique( vec.atom.type.mutation ),
                                                                  unique( vec.atom.type.non.mutation ) ) )
    df.atom.type.mutation.nonmutation[ , 1 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 1 ] )
    df.atom.type.mutation.nonmutation[ , 2 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 2 ] )
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- NULL
    for( j in 1:nrow( df.atom.type.mutation.nonmutation ) )
    {
      r1 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 1 ], ]$rvdw
      r2 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 2 ], ]$rvdw
      df.atom.type.mutation.nonmutation.r1.r2.0.9 <- rbind( df.atom.type.mutation.nonmutation.r1.r2.0.9,
                                                            c( df.atom.type.mutation.nonmutation[ j, ],
                                                               r1 + r2 - overlap.cutoff ) )
    }
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- data.frame( df.atom.type.mutation.nonmutation.r1.r2.0.9 )
    colnames( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) <- c( "mutation",
                                                                  "nonmutation",
                                                                  "value" )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$value <- as.numeric( as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$value ) ) )
    
    mat.r1.r2.0.9 <- matrix( 0,
                             length( vec.atom.type.mutation ),
                             length( vec.atom.type.non.mutation ) )
    
    for( j in 1:nrow( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) )
    {
      vec.i <- which( vec.atom.type.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$mutation )
      vec.j <- which( vec.atom.type.non.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$nonmutation )
      
      mat.r1.r2.0.9[ vec.i, vec.j ] <- df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$value
    }
    
    # Get number of clashes
    mat.clash <- mat.dist < mat.r1.r2.0.9
    
    # Return output
    # vec.prefusion.no.clashes <- c( vec.prefusion.no.clashes,
    #                                sum( mat.clash ) )
    
    sum( mat.clash )
  }
  # vec.prefusion.no.clashes <- unlist( vec.prefusion.no.clashes )
  
  df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes <- cbind( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain,
                                                                               vec.prefusion.no.clashes )
  colnames( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes ) <- c( colnames( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ),
                                                                                       "no_clash_prefusion" )
  write.table( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes,
               file = paste( output_prefix,
                             "_prefusion.txt",
                             sep = "" ) )
  
  #==================================================
  #==================================================
  #==================================================
  # Postfusion conformation
  
  #==================================================
  # Load original PDB and output altered one removing all non-ATOM entries
  
  # Load input PDB
  pdb.postfusion.raw <- read.pdb( file.pdb.postfusion )
  pdb.postfusion.raw.atom <- pdb.postfusion.raw$atom
  
  # Save new altered PDB
  pdb.postfusion.raw.atom.altered <- pdb.postfusion.raw.atom[ pdb.postfusion.raw.atom$type == "ATOM", ]
  file.pdb.postfusion.altered <- paste( output_prefix,
                                        "_postfusion.pdb",
                                        sep = "" )
  write.pdb( file = file.pdb.postfusion.altered,
             xyz = as.vector( t( cbind( pdb.postfusion.raw.atom.altered$x,
                                        pdb.postfusion.raw.atom.altered$y,
                                        pdb.postfusion.raw.atom.altered$z ) ) ),
             type = pdb.postfusion.raw.atom.altered$type,
             eleno = pdb.postfusion.raw.atom.altered$eleno,
             elety = pdb.postfusion.raw.atom.altered$elety,
             alt = pdb.postfusion.raw.atom.altered$alt,
             resid = pdb.postfusion.raw.atom.altered$resid,
             chain = pdb.postfusion.raw.atom.altered$chain,
             resno = pdb.postfusion.raw.atom.altered$resno,
             insert = pdb.postfusion.raw.atom.altered$insert,
             o = pdb.postfusion.raw.atom.altered$o,
             b = pdb.postfusion.raw.atom.altered$b )
  
  pdb.postfusion <- read.pdb( file.pdb.postfusion.altered )
  pdb.postfusion.atom <- pdb.postfusion$atom
  
  #==================================================
  # Generate data frame that saves chain and resno
  df.chain.resno.resid.postfusion <- NULL
  for( j in 1:nrow( pdb.postfusion.atom ) )
  {
    var.chain <- pdb.postfusion.atom[ j, ]$chain
    if( is.na( pdb.postfusion.atom[ j, ]$insert ) )
    {
      var.resno <- as.character( pdb.postfusion.atom[ j, ]$resno )
    } else
    {
      var.resno <- paste( pdb.postfusion.atom[ j, ]$resno,
                          pdb.postfusion.atom[ j, ]$insert,
                          sep = "" )
    }
    var.resid <- pdb.postfusion.atom[ j, ]$resid
    
    if( length( df.chain.resno.resid.postfusion[ df.chain.resno.resid.postfusion[ , 1 ] == var.chain &
                                                 df.chain.resno.resid.postfusion[ , 2 ] == var.resno, ] ) == 0 )
    {
      df.chain.resno.resid.postfusion <- rbind( df.chain.resno.resid.postfusion,
                                                c( var.chain,
                                                   var.resno,
                                                   var.resid ) )
    }
  }
  df.chain.resno.resid.postfusion <- data.frame( df.chain.resno.resid.postfusion )
  colnames( df.chain.resno.resid.postfusion ) <- c( "chain_postfusion",
                                                    "resno_postfusion",
                                                    "resid_postfusion" )
  df.chain.resno.resid.postfusion$chain_postfusion <- as.vector( df.chain.resno.resid.postfusion$chain_postfusion )
  df.chain.resno.resid.postfusion$resno_postfusion <- as.vector( df.chain.resno.resid.postfusion$resno_postfusion )
  df.chain.resno.resid.postfusion$resid_postfusion <- as.vector( aa321( df.chain.resno.resid.postfusion$resid_postfusion ) )
  
  #==================================================
  # Calculate torsion angles for PDB using DSSP
  
  # Determine torsion angles (and secondary structure) using DSSP
  list.ss.dssp.postfusion <- dssp( pdb.postfusion,
                                   resno = TRUE )
  vec.phi.postfusion <- list.ss.dssp.postfusion$phi
  vec.phi.postfusion[ which( vec.phi.postfusion == 360 ) ] <- NA
  
  vec.psi.postfusion <- list.ss.dssp.postfusion$psi
  vec.psi.postfusion[ which( vec.psi.postfusion == 360 ) ] <- NA
  
  # Select torsion angles of interest
  df.postfusion.phi.psi <- data.frame( cbind( vec.phi.postfusion,
                                              vec.psi.postfusion ) )
  
  #==================================================
  # Combine PDB infrmation with torsion angles
  df.chain.resno.resid.postfusion.phi.psi <- data.frame( cbind( df.chain.resno.resid.postfusion,
                                                                df.postfusion.phi.psi ) )
  colnames( df.chain.resno.resid.postfusion.phi.psi ) <- c( colnames( df.chain.resno.resid.postfusion ),
                                                            "phi_postfusion",
                                                            "psi_postfusion" )
  
  #==================================================
  # Identify if angles are in trans proline conformation
  vec.phi.psi.postfusion.trans.proline <- NULL
  for( i in 1:nrow( df.chain.resno.resid.postfusion.phi.psi) )
  {
    phi <- df.chain.resno.resid.postfusion.phi.psi[ i, ]$phi_postfusion
    psi <- df.chain.resno.resid.postfusion.phi.psi[ i, ]$psi_postfusion
    
    if( is.na( phi ) | is.na( psi ) )
    {
      vec.phi.psi.postfusion.trans.proline <- c( vec.phi.psi.postfusion.trans.proline,
                                                 NA )
    } else
    {
      phi.tmp <- abs( phi - vec.potential.x )
      phi.tmp.choice <- which( phi.tmp == min( phi.tmp ) )
      if( length( phi.tmp.choice ) > 1 )
      {
        # phi.tmp.choice <- sample( phi.tmp.choice,
        #                           1 )
        phi.tmp.choice <- min( phi.tmp.choice )
      }
      phi.round <- vec.potential.x[ phi.tmp.choice ]
      phi.round
      
      psi.tmp <- abs( psi - vec.potential.x )
      psi.tmp.choice <- which( psi.tmp == min( psi.tmp ) )
      if( length( psi.tmp.choice ) > 1 )
      {
        # psi.tmp.choice <- sample( psi.tmp.choice,
        #                           1 )
        psi.tmp.choice <- min( psi.tmp.choice )
      }
      psi.round <- vec.potential.x[ psi.tmp.choice ]
      psi.round
      
      if( nrow( df.x.y[ df.x.y$x == phi.round &
                        df.x.y$y == psi.round, ] ) > 0 )
      {
        vec.phi.psi.postfusion.trans.proline <- c( vec.phi.psi.postfusion.trans.proline,
                                                   1 )
      } else
      {
        vec.phi.psi.postfusion.trans.proline <- c( vec.phi.psi.postfusion.trans.proline,
                                                   0 )
      }
    }
  }
  
  # Add column to data frame
  df.chain.resno.resid.postfusion.phi.psi.angle <- cbind( df.chain.resno.resid.postfusion.phi.psi,
                                                          vec.phi.psi.postfusion.trans.proline )
  colnames( df.chain.resno.resid.postfusion.phi.psi.angle ) <- c( colnames( df.chain.resno.resid.postfusion.phi.psi ),
                                                                  "angle_trans_proline_postfusion" )
  
  #==================================================
  # Secondary structure
  
  # Add extra column to add helix conformations
  df.chain.resno.resid.postfusion.phi.psi.angle.helix <- cbind( df.chain.resno.resid.postfusion.phi.psi.angle,
                                                                rep( 0,
                                                                     nrow( df.chain.resno.resid.postfusion.phi.psi.angle ) ),
                                                                rep( 0,
                                                                     nrow( df.chain.resno.resid.postfusion.phi.psi.angle ) ),
                                                                rep( 0,
                                                                     nrow( df.chain.resno.resid.postfusion.phi.psi.angle ) ) )
  colnames( df.chain.resno.resid.postfusion.phi.psi.angle.helix ) <- c( colnames( df.chain.resno.resid.postfusion.phi.psi.angle ),
                                                                        "helix_postfusion",
                                                                        "helix_termini_postfusion",
                                                                        "DSSP_postfusion" )
  
  
  vec.ind.subset.chain <- as.numeric( which( list.ss.dssp.postfusion$helix$chain == postfusion.chain ) )
  
  vec.helix.start <- list.ss.dssp.postfusion$helix$start[ vec.ind.subset.chain ]
  vec.helix.end <- list.ss.dssp.postfusion$helix$end[ vec.ind.subset.chain ]
  vec.index.helix <- NULL
  if( length( vec.helix.start ) > 0 )
  {
    for( i in 1:length( vec.helix.start ) )
    {
      vec.index.helix <- c( vec.index.helix,
                            vec.helix.start[ i ]:vec.helix.end[ i ] )
    }
    
    vec.resno.helix <- df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% as.character( vec.index.helix ), ]$resno_postfusion
    # Add to data frame
    df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% vec.resno.helix, ]$helix_postfusion <- 1
    # df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% c( vec.helix.start, vec.helix.end ), ]$helix_termini_postfusion <- 1
    df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% as.character( sort( as.vector( c( vec.helix.start, vec.helix.end ) ) ) ), ]$helix_termini_postfusion <- 1
  }
  
  
  # Add DSSP secondary structure to data frame
  vec.dssp <- list.ss.dssp.postfusion$sse
  vec.dssp[ vec.dssp == " " ] = "L"
  df.chain.resno.resid.postfusion.phi.psi.angle.helix$DSSP_postfusion <- vec.dssp
  
  
  #==================================================
  # Run pymol to mutate all positions
  print( "Run pymol to mutate all positions in the postfusion conformation" )
  vec.postfusion.no.clashes <- NULL
  df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain <- df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$chain_postfusion == postfusion.chain, ]
  # for( i in 1:nrow( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ) )
  vec.postfusion.no.clashes <- foreach( i = 1:nrow( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ), .combine = c ) %dopar%
  {
    # start.time <- proc.time()
    # print( paste( i,
    #               nrow( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ),
    #               sep = "/" ) )
    
    chain <- df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain[ i, 1 ]
    resno <- df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain[ i, 2 ]
    selection <- paste( chain,
                        "/",
                        resno,
                        "/",
                        sep = "" )
    
    file.pdb.postfusion.mutation <- paste( unlist( strsplit( file.pdb.postfusion,
                                                             "\\.pdb" ) ),
                                           "_",
                                           chain,
                                           resno,
                                           ".pdb",
                                           sep = "" )
    
    #==================================================
    # Generate mutated PDB file
    command <- paste( "pymol",
                      "-qc",
                      # "~/Documents/work/vrc/sbis_projects/proline_predictor/scripts/mutate.py",
                      # "mutate.py",
                      file.exe.mutate,
                      "--",
                      file.pdb.postfusion,
                      file.pdb.postfusion.mutation,
                      selection,
                      "PRO",
                      distance_angstroem )
    # Run python (PyMol) command
    pymol.out <- system( command,
                         intern = TRUE )
    if( sum( grepl( "Error",
                    pymol.out ) ) > 0 )
    {
      # vec.postfusion.no.clashes <- c( vec.postfusion.no.clashes,
      #                                 0 )
      # next
      
      return( 0 )
    }
    
    #==================================================
    # Analyze mutated protein structure
    
    # Load mutated structure
    pdb.mutant <- read.pdb( file.pdb.postfusion.mutation )
    pdb.mutant.atom <- pdb.mutant$atom
    
    # Select mutated residue
    # Do we have an insertion?
    if( is.letter( resno ) )
    {
      resno.new <- NULL
      insert.new <- NULL
      for( char in unlist( strsplit( resno,
                                     "" ) ) )
      {
        if( is.letter( char ) )
        {
          insert.new <- c( insert.new,
                           char )
        } else
        {
          resno.new <- c( resno.new,
                          char )
        }
      }
      
      resno.new <- as.numeric( paste( resno.new,
                                      collapse = "" ) )
      insert.new <- paste( insert.new,
                           collapse = "" )
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno.new &
                                                     pdb.mutant.atom$insert == insert.new &
                                                     !( is.na( pdb.mutant.atom$insert ) ), ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno.new &
                                                            pdb.mutant.atom$insert == insert.new &
                                                            !( is.na( pdb.mutant.atom$insert ) ) ), ]
      
    } else
    {
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno, ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno ), ]
    }
    # Select xyz coordinates for non-backbone atoms of mutated residue 
    pdb.mutant.atom.mutation.no.backbone <- pdb.mutant.atom.mutation[ which( !( pdb.mutant.atom.mutation$elety %in% c( "N",
                                                                                                                       "CA",
                                                                                                                       "C",
                                                                                                                       "O" ) ) ), ]
    pdb.mutant.atom.mutation.no.backbone.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.mutation.no.backbone[ , 9:11 ] ) ) )
    
    # Select xyz coordinates for all other residues
    pdb.mutant.atom.non.mutation.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.non.mutation[ , 9:11 ] ) ) )
    
    # Calculate distance matrix
    mat.dist <- dist.xyz( pdb.mutant.atom.mutation.no.backbone.xyz,
                          pdb.mutant.atom.non.mutation.xyz )
    
    # Generate vectors representing atom types for mutated residue and rest
    vec.atom.type.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.mutation.no.backbone$elety ) )
    {
      elety <- pdb.mutant.atom.mutation.no.backbone$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.mutation <- c( vec.atom.type.mutation,
                                   elety.split[ which( unlist( lapply( elety.split,
                                                                       is.letter ) ) )[ 1 ] ] )
    }
    
    vec.atom.type.non.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.non.mutation$elety ) )
    {
      elety <- pdb.mutant.atom.non.mutation$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.non.mutation <- c( vec.atom.type.non.mutation,
                                       elety.split[ which( unlist( lapply( elety.split,
                                                                           is.letter ) ) )[ 1 ] ] )
    }
    
    # Calculate MolProbity clash constraint
    df.atom.type.mutation.nonmutation <- data.frame( expand.grid( unique( vec.atom.type.mutation ),
                                                                  unique( vec.atom.type.non.mutation ) ) )
    df.atom.type.mutation.nonmutation[ , 1 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 1 ] )
    df.atom.type.mutation.nonmutation[ , 2 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 2 ] )
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- NULL
    for( j in 1:nrow( df.atom.type.mutation.nonmutation ) )
    {
      r1 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 1 ], ]$rvdw
      r2 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 2 ], ]$rvdw
      df.atom.type.mutation.nonmutation.r1.r2.0.9 <- rbind( df.atom.type.mutation.nonmutation.r1.r2.0.9,
                                                            c( df.atom.type.mutation.nonmutation[ j, ],
                                                               r1 + r2 - overlap.cutoff ) )
    }
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- data.frame( df.atom.type.mutation.nonmutation.r1.r2.0.9 )
    colnames( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) <- c( "mutation",
                                                                  "nonmutation",
                                                                  "value" )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$value <- as.numeric( as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$value ) ) )
    
    mat.r1.r2.0.9 <- matrix( 0,
                             length( vec.atom.type.mutation ),
                             length( vec.atom.type.non.mutation ) )
    
    for( j in 1:nrow( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) )
    {
      vec.i <- which( vec.atom.type.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$mutation )
      vec.j <- which( vec.atom.type.non.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$nonmutation )
      
      mat.r1.r2.0.9[ vec.i, vec.j ] <- df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$value
    }
    
    # Get number of clashes
    mat.clash <- mat.dist < mat.r1.r2.0.9
    
    # vec.postfusion.no.clashes <- c( vec.postfusion.no.clashes,
    #                                 sum( mat.clash ) )
    # print( vec.postfusion.no.clashes )
    # print( ( proc.time() - start.time )[ 3 ] )
    sum( mat.clash )
  }
  # vec.postfusion.no.clashes <- unlist( vec.postfusion.no.clashes )
  
  df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes <- cbind( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain,
                                                                                vec.postfusion.no.clashes )
  colnames( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes ) <- c( colnames( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ),
                                                                                        "no_clash_postfusion" )
  write.table( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes,
               file = paste( output_prefix,
                             "_postfusion.txt",
                             sep = "" ) )
} else if( option == 2 )
{
  #==================================================
  #==================================================
  #==================================================
  # Prefusion conformation
  
  #==================================================
  # Load original PDB and output altered one removing all non-ATOM entries
  
  # Load input PDB
  pdb.prefusion.raw <- read.pdb( file.pdb.prefusion )
  pdb.prefusion.raw.atom <- pdb.prefusion.raw$atom
  
  # Save new altered PDB
  pdb.prefusion.raw.atom.altered <- pdb.prefusion.raw.atom[ pdb.prefusion.raw.atom$type == "ATOM", ]
  file.pdb.prefusion.altered <- paste( output_prefix,
                                       "_prefusion.pdb",
                                       sep = "" )
  write.pdb( file = file.pdb.prefusion.altered,
             xyz = as.vector( t( cbind( pdb.prefusion.raw.atom.altered$x,
                                        pdb.prefusion.raw.atom.altered$y,
                                        pdb.prefusion.raw.atom.altered$z ) ) ),
             type = pdb.prefusion.raw.atom.altered$type,
             eleno = pdb.prefusion.raw.atom.altered$eleno,
             elety = pdb.prefusion.raw.atom.altered$elety,
             alt = pdb.prefusion.raw.atom.altered$alt,
             resid = pdb.prefusion.raw.atom.altered$resid,
             chain = pdb.prefusion.raw.atom.altered$chain,
             resno = pdb.prefusion.raw.atom.altered$resno,
             insert = pdb.prefusion.raw.atom.altered$insert,
             o = pdb.prefusion.raw.atom.altered$o,
             b = pdb.prefusion.raw.atom.altered$b )
  
  pdb.prefusion <- read.pdb( file.pdb.prefusion.altered )
  pdb.prefusion.atom <- pdb.prefusion$atom
  
  #==================================================
  # Generate data frame that saves chain and resno
  df.chain.resno.resid.prefusion <- NULL
  for( j in 1:nrow( pdb.prefusion.atom ) )
  {
    var.chain <- pdb.prefusion.atom[ j, ]$chain
    if( is.na( pdb.prefusion.atom[ j, ]$insert ) )
    {
      var.resno <- as.character( pdb.prefusion.atom[ j, ]$resno )
    } else
    {
      var.resno <- paste( pdb.prefusion.atom[ j, ]$resno,
                          pdb.prefusion.atom[ j, ]$insert,
                          sep = "" )
    }
    var.resid <- pdb.prefusion.atom[ j, ]$resid
    
    if( length( df.chain.resno.resid.prefusion[ df.chain.resno.resid.prefusion[ , 1 ] == var.chain &
                                                df.chain.resno.resid.prefusion[ , 2 ] == var.resno, ] ) == 0 )
    {
      df.chain.resno.resid.prefusion <- rbind( df.chain.resno.resid.prefusion,
                                               c( var.chain,
                                                  var.resno,
                                                  var.resid ) )
    }
  }
  df.chain.resno.resid.prefusion <- data.frame( df.chain.resno.resid.prefusion )
  colnames( df.chain.resno.resid.prefusion ) <- c( "chain_prefusion",
                                                   "resno_prefusion",
                                                   "resid_prefusion" )
  df.chain.resno.resid.prefusion$chain_prefusion <- as.vector( df.chain.resno.resid.prefusion$chain_prefusion )
  df.chain.resno.resid.prefusion$resno_prefusion <- as.vector( df.chain.resno.resid.prefusion$resno_prefusion )
  df.chain.resno.resid.prefusion$resid_prefusion <- as.vector( aa321( df.chain.resno.resid.prefusion$resid_prefusion ) )
  
  #==================================================
  # Calculate torsion angles for PDB using DSSP
  
  # Determine torsion angles (and secondary structure) using DSSP
  list.ss.dssp.prefusion <- dssp( pdb.prefusion,
                                  resno = TRUE ) 
  vec.phi.prefusion <- list.ss.dssp.prefusion$phi
  vec.phi.prefusion[ which( vec.phi.prefusion == 360 ) ] <- NA
  
  vec.psi.prefusion <- list.ss.dssp.prefusion$psi
  vec.psi.prefusion[ which( vec.psi.prefusion == 360 ) ] <- NA
  
  # Select torsion angles of interest
  df.prefusion.phi.psi <- data.frame( cbind( vec.phi.prefusion,
                                             vec.psi.prefusion ) )
  
  #==================================================
  # Combine PDB information with torsion angles
  df.chain.resno.resid.prefusion.phi.psi <- data.frame( cbind( df.chain.resno.resid.prefusion,
                                                               df.prefusion.phi.psi ) )
  colnames( df.chain.resno.resid.prefusion.phi.psi ) <- c( colnames( df.chain.resno.resid.prefusion ),
                                                           "phi_prefusion",
                                                           "psi_prefusion" )
  
  #==================================================
  # Identify if angles are in trans proline conformation
  vec.phi.psi.prefusion.trans.proline <- NULL
  for( i in 1:nrow( df.chain.resno.resid.prefusion.phi.psi) )
  {
    phi <- df.chain.resno.resid.prefusion.phi.psi[ i, ]$phi_prefusion
    psi <- df.chain.resno.resid.prefusion.phi.psi[ i, ]$psi_prefusion
    
    if( is.na( phi ) | is.na( psi ) )
    {
      vec.phi.psi.prefusion.trans.proline <- c( vec.phi.psi.prefusion.trans.proline,
                                                NA )
    } else
    {
      phi.tmp <- abs( phi - vec.potential.x )
      phi.tmp.choice <- which( phi.tmp == min( phi.tmp ) )
      if( length( phi.tmp.choice ) > 1 )
      {
        # phi.tmp.choice <- sample( phi.tmp.choice,
        #                           1 )
        phi.tmp.choice <- min( phi.tmp.choice )
      }
      phi.round <- vec.potential.x[ phi.tmp.choice ]
      phi.round
      
      psi.tmp <- abs( psi - vec.potential.x )
      psi.tmp.choice <- which( psi.tmp == min( psi.tmp ) )
      if( length( psi.tmp.choice ) > 1 )
      {
        # psi.tmp.choice <- sample( psi.tmp.choice,
        #                           1 )
        psi.tmp.choice <- min( psi.tmp.choice )
      }
      psi.round <- vec.potential.x[ psi.tmp.choice ]
      psi.round
      
      if( nrow( df.x.y[ df.x.y$x == phi.round &
                        df.x.y$y == psi.round, ] ) > 0 )
      {
        vec.phi.psi.prefusion.trans.proline <- c( vec.phi.psi.prefusion.trans.proline,
                                                  1 )
      } else
      {
        vec.phi.psi.prefusion.trans.proline <- c( vec.phi.psi.prefusion.trans.proline,
                                                  0 )
      }
    }
  }
  
  # Add column to data frame
  df.chain.resno.resid.prefusion.phi.psi.angle <- cbind( df.chain.resno.resid.prefusion.phi.psi,
                                                         vec.phi.psi.prefusion.trans.proline )
  colnames( df.chain.resno.resid.prefusion.phi.psi.angle ) <- c( colnames( df.chain.resno.resid.prefusion.phi.psi ),
                                                                 "angle_trans_proline_prefusion" )
  
  
  
  #==================================================
  # Secondary structure and PHI and PSI angles
  
  # Add extra column to add helix conformations
  df.chain.resno.resid.prefusion.phi.psi.angle.helix <- cbind( df.chain.resno.resid.prefusion.phi.psi.angle,
                                                               rep( 0,
                                                                    nrow( df.chain.resno.resid.prefusion.phi.psi.angle ) ),
                                                               rep( 0,
                                                                    nrow( df.chain.resno.resid.prefusion.phi.psi.angle ) ),
                                                               rep( 0,
                                                                    nrow( df.chain.resno.resid.prefusion.phi.psi.angle ) ) )
  colnames( df.chain.resno.resid.prefusion.phi.psi.angle.helix ) <- c( colnames( df.chain.resno.resid.prefusion.phi.psi.angle ),
                                                                       "helix_prefusion",
                                                                       "helix_termini_prefusion",
                                                                       "DSSP_prefusion" )
  
  
  # Select subsets where chain of interest if given
  vec.ind.subset.chain <- as.numeric( which( list.ss.dssp.prefusion$helix$chain == prefusion.chain ) )
  
  vec.helix.start <- list.ss.dssp.prefusion$helix$start[ vec.ind.subset.chain ]
  vec.helix.end <- list.ss.dssp.prefusion$helix$end[ vec.ind.subset.chain ]
  vec.index.helix <- NULL
  if( length( vec.helix.start ) > 0 )
  {
    for( i in 1:length( vec.helix.start ) )
    {
      vec.index.helix <- c( vec.index.helix,
                            vec.helix.start[ i ]:vec.helix.end[ i ] )
    }
    
    vec.resno.helix <- df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% as.character( vec.index.helix ), ]$resno_prefusion
    
    #==================================================
    # Add to data frame
    df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% vec.resno.helix, ]$helix_prefusion <- 1
    # df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% c( vec.helix.start, vec.helix.end ), ]$helix_termini_prefusion <- 1
    df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$resno_prefusion %in% as.character( sort( as.vector( c( vec.helix.start, vec.helix.end ) ) ) ), ]$helix_termini_prefusion <- 1
  }
  
  # Add DSSP secondary structure to data frame
  vec.dssp <- list.ss.dssp.prefusion$sse
  vec.dssp[ vec.dssp == " " ] = "L"
  df.chain.resno.resid.prefusion.phi.psi.angle.helix$DSSP_prefusion <- vec.dssp
  
  #==================================================
  # Run pymol to mutate all positions
  print( "Run pymol to mutate all positions in the prefusion conformation" )
  df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain <- df.chain.resno.resid.prefusion.phi.psi.angle.helix[ df.chain.resno.resid.prefusion.phi.psi.angle.helix$chain_prefusion == prefusion.chain, ]
  
  # vec.prefusion.no.clashes <- NULL
  # for( i in 1:nrow( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ) )
  # {
  #   start.time <- proc.time()
  #   print( paste( i,
  #                 nrow( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ),
  #                 sep = "/" ) )
  vec.prefusion.no.clashes <- foreach( i = 1:nrow( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ), .combine = c ) %dopar%
  {
    chain <- df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain[ i, 1 ]
    resno <- df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain[ i, 2 ]
    selection <- paste( chain,
                        "/",
                        resno,
                        "/",
                        sep = "" )
    
    file.pdb.prefusion.mutation <- paste( unlist( strsplit( file.pdb.prefusion,
                                                            "\\.pdb" ) ),
                                          "_",
                                          chain,
                                          resno,
                                          ".pdb",
                                          sep = "" )
    
    #==================================================
    # Generate mutated PDB file
    command <- paste( "pymol",
                      "-qc",
                      # "~/Documents/work/vrc/sbis_projects/proline_predictor/scripts/mutate.py",
                      # "mutate.py",
                      file.exe.mutate,
                      "--",
                      file.pdb.prefusion,
                      file.pdb.prefusion.mutation,
                      selection,
                      "PRO",
                      distance_angstroem )
    # Run python (PyMol) command
    system( command,
            intern = TRUE )
    
    #==================================================
    # Analyze mutated protein structure
    
    # Load mutated structure
    pdb.mutant <- read.pdb( file.pdb.prefusion.mutation )
    pdb.mutant.atom <- pdb.mutant$atom
    
    # Select mutated residue
    # Do we have an insertion?
    if( is.letter( resno ) )
    {
      resno.new <- NULL
      insert.new <- NULL
      for( char in unlist( strsplit( resno,
                                     "" ) ) )
      {
        if( is.letter( char ) )
        {
          insert.new <- c( insert.new,
                           char )
        } else
        {
          resno.new <- c( resno.new,
                          char )
        }
      }
      
      resno.new <- as.numeric( paste( resno.new,
                                      collapse = "" ) )
      insert.new <- paste( insert.new,
                           collapse = "" )
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno.new &
                                                     pdb.mutant.atom$insert == insert.new &
                                                     !( is.na( pdb.mutant.atom$insert ) ), ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno.new &
                                                            pdb.mutant.atom$insert == insert.new &
                                                            !( is.na( pdb.mutant.atom$insert ) ) ), ]
      
    } else
    {
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno, ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno ), ]
    }
    # Select xyz coordinates for non-backbone atoms of mutated residue
    pdb.mutant.atom.mutation.no.backbone <- pdb.mutant.atom.mutation[ which( !( pdb.mutant.atom.mutation$elety %in% c( "N",
                                                                                                                       "CA",
                                                                                                                       "C",
                                                                                                                       "O" ) ) ), ]
    pdb.mutant.atom.mutation.no.backbone.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.mutation.no.backbone[ , 9:11 ] ) ) )
    
    # Select xyz coordinates for all other residues
    pdb.mutant.atom.non.mutation.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.non.mutation[ , 9:11 ] ) ) )
    
    # Calculate distance matrix
    mat.dist <- dist.xyz( pdb.mutant.atom.mutation.no.backbone.xyz,
                          pdb.mutant.atom.non.mutation.xyz )
    
    # Generate vectors representing atom types for mutated residue and rest
    vec.atom.type.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.mutation.no.backbone$elety ) )
    {
      elety <- pdb.mutant.atom.mutation.no.backbone$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.mutation <- c( vec.atom.type.mutation,
                                   elety.split[ which( unlist( lapply( elety.split,
                                                                       is.letter ) ) )[ 1 ] ] )
    }
    
    vec.atom.type.non.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.non.mutation$elety ) )
    {
      elety <- pdb.mutant.atom.non.mutation$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.non.mutation <- c( vec.atom.type.non.mutation,
                                       elety.split[ which( unlist( lapply( elety.split,
                                                                           is.letter ) ) )[ 1 ] ] )
    }
    
    # Calculate MolProbity clash constraint
    df.atom.type.mutation.nonmutation <- data.frame( expand.grid( unique( vec.atom.type.mutation ),
                                                                  unique( vec.atom.type.non.mutation ) ) )
    df.atom.type.mutation.nonmutation[ , 1 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 1 ] )
    df.atom.type.mutation.nonmutation[ , 2 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 2 ] )
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- NULL
    for( j in 1:nrow( df.atom.type.mutation.nonmutation ) )
    {
      r1 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 1 ], ]$rvdw
      r2 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 2 ], ]$rvdw
      df.atom.type.mutation.nonmutation.r1.r2.0.9 <- rbind( df.atom.type.mutation.nonmutation.r1.r2.0.9,
                                                            c( df.atom.type.mutation.nonmutation[ j, ],
                                                               r1 + r2 - overlap.cutoff ) )
    }
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- data.frame( df.atom.type.mutation.nonmutation.r1.r2.0.9 )
    colnames( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) <- c( "mutation",
                                                                  "nonmutation",
                                                                  "value" )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$value <- as.numeric( as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$value ) ) )
    
    mat.r1.r2.0.9 <- matrix( 0,
                             length( vec.atom.type.mutation ),
                             length( vec.atom.type.non.mutation ) )
    
    for( j in 1:nrow( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) )
    {
      vec.i <- which( vec.atom.type.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$mutation )
      vec.j <- which( vec.atom.type.non.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$nonmutation )
      
      mat.r1.r2.0.9[ vec.i, vec.j ] <- df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$value
    }
    
    # Get number of clashes
    mat.clash <- mat.dist < mat.r1.r2.0.9
    
    # Return output
    # vec.prefusion.no.clashes <- c( vec.prefusion.no.clashes,
    #                                sum( mat.clash ) )
    
    sum( mat.clash )
  }
  # vec.prefusion.no.clashes <- unlist( vec.prefusion.no.clashes )
  
  df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes <- cbind( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain,
                                                                               vec.prefusion.no.clashes )
  colnames( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes ) <- c( colnames( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain ),
                                                                                       "no_clash_prefusion" )
  write.table( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes,
               file = paste( output_prefix,
                             "_prefusion.txt",
                             sep = "" ) )
} else if( option == 3 )
{
  print( "Option 3" )
  #==================================================
  # Load original PDB and output altered one removing all non-ATOM entries
  
  # Load input PDB
  pdb.postfusion.raw <- read.pdb( file.pdb.postfusion )
  pdb.postfusion.raw.atom <- pdb.postfusion.raw$atom
  
  # Save new altered PDB
  pdb.postfusion.raw.atom.altered <- pdb.postfusion.raw.atom[ pdb.postfusion.raw.atom$type == "ATOM", ]
  file.pdb.postfusion.altered <- paste( output_prefix,
                                        "_postfusion.pdb",
                                        sep = "" )
  write.pdb( file = file.pdb.postfusion.altered,
             xyz = as.vector( t( cbind( pdb.postfusion.raw.atom.altered$x,
                                        pdb.postfusion.raw.atom.altered$y,
                                        pdb.postfusion.raw.atom.altered$z ) ) ),
             type = pdb.postfusion.raw.atom.altered$type,
             eleno = pdb.postfusion.raw.atom.altered$eleno,
             elety = pdb.postfusion.raw.atom.altered$elety,
             alt = pdb.postfusion.raw.atom.altered$alt,
             resid = pdb.postfusion.raw.atom.altered$resid,
             chain = pdb.postfusion.raw.atom.altered$chain,
             resno = pdb.postfusion.raw.atom.altered$resno,
             insert = pdb.postfusion.raw.atom.altered$insert,
             o = pdb.postfusion.raw.atom.altered$o,
             b = pdb.postfusion.raw.atom.altered$b )
  
  pdb.postfusion <- read.pdb( file.pdb.postfusion.altered )
  pdb.postfusion.atom <- pdb.postfusion$atom
  
  #==================================================
  # Generate data frame that saves chain and resno
  df.chain.resno.resid.postfusion <- NULL
  for( j in 1:nrow( pdb.postfusion.atom ) )
  {
    var.chain <- pdb.postfusion.atom[ j, ]$chain
    if( is.na( pdb.postfusion.atom[ j, ]$insert ) )
    {
      var.resno <- as.character( pdb.postfusion.atom[ j, ]$resno )
    } else
    {
      var.resno <- paste( pdb.postfusion.atom[ j, ]$resno,
                          pdb.postfusion.atom[ j, ]$insert,
                          sep = "" )
    }
    var.resid <- pdb.postfusion.atom[ j, ]$resid
    
    if( length( df.chain.resno.resid.postfusion[ df.chain.resno.resid.postfusion[ , 1 ] == var.chain &
                                                 df.chain.resno.resid.postfusion[ , 2 ] == var.resno, ] ) == 0 )
    {
      df.chain.resno.resid.postfusion <- rbind( df.chain.resno.resid.postfusion,
                                                c( var.chain,
                                                   var.resno,
                                                   var.resid ) )
    }
  }
  df.chain.resno.resid.postfusion <- data.frame( df.chain.resno.resid.postfusion )
  colnames( df.chain.resno.resid.postfusion ) <- c( "chain_postfusion",
                                                    "resno_postfusion",
                                                    "resid_postfusion" )
  df.chain.resno.resid.postfusion$chain_postfusion <- as.vector( df.chain.resno.resid.postfusion$chain_postfusion )
  df.chain.resno.resid.postfusion$resno_postfusion <- as.vector( df.chain.resno.resid.postfusion$resno_postfusion )
  df.chain.resno.resid.postfusion$resid_postfusion <- as.vector( aa321( df.chain.resno.resid.postfusion$resid_postfusion ) )
  
  #==================================================
  # Calculate torsion angles for PDB using DSSP
  
  # Determine torsion angles (and secondary structure) using DSSP
  list.ss.dssp.postfusion <- dssp( pdb.postfusion,
                                   resno = TRUE )
  vec.phi.postfusion <- list.ss.dssp.postfusion$phi
  vec.phi.postfusion[ which( vec.phi.postfusion == 360 ) ] <- NA
  
  vec.psi.postfusion <- list.ss.dssp.postfusion$psi
  vec.psi.postfusion[ which( vec.psi.postfusion == 360 ) ] <- NA
  
  # Select torsion angles of interest
  df.postfusion.phi.psi <- data.frame( cbind( vec.phi.postfusion,
                                              vec.psi.postfusion ) )
  
  #==================================================
  # Combine PDB infrmation with torsion angles
  df.chain.resno.resid.postfusion.phi.psi <- data.frame( cbind( df.chain.resno.resid.postfusion,
                                                                df.postfusion.phi.psi ) )
  colnames( df.chain.resno.resid.postfusion.phi.psi ) <- c( colnames( df.chain.resno.resid.postfusion ),
                                                            "phi_postfusion",
                                                            "psi_postfusion" )
  
  #==================================================
  # Identify if angles are in trans proline conformation
  vec.phi.psi.postfusion.trans.proline <- NULL
  for( i in 1:nrow( df.chain.resno.resid.postfusion.phi.psi) )
  {
    phi <- df.chain.resno.resid.postfusion.phi.psi[ i, ]$phi_postfusion
    psi <- df.chain.resno.resid.postfusion.phi.psi[ i, ]$psi_postfusion
    
    if( is.na( phi ) | is.na( psi ) )
    {
      vec.phi.psi.postfusion.trans.proline <- c( vec.phi.psi.postfusion.trans.proline,
                                                 NA )
    } else
    {
      phi.tmp <- abs( phi - vec.potential.x )
      phi.tmp.choice <- which( phi.tmp == min( phi.tmp ) )
      if( length( phi.tmp.choice ) > 1 )
      {
        # phi.tmp.choice <- sample( phi.tmp.choice,
        #                           1 )
        phi.tmp.choice <- min( phi.tmp.choice )
      }
      phi.round <- vec.potential.x[ phi.tmp.choice ]
      phi.round
      
      psi.tmp <- abs( psi - vec.potential.x )
      psi.tmp.choice <- which( psi.tmp == min( psi.tmp ) )
      if( length( psi.tmp.choice ) > 1 )
      {
        # psi.tmp.choice <- sample( psi.tmp.choice,
        #                           1 )
        psi.tmp.choice <- min( psi.tmp.choice )
      }
      psi.round <- vec.potential.x[ psi.tmp.choice ]
      psi.round
      
      if( nrow( df.x.y[ df.x.y$x == phi.round &
                        df.x.y$y == psi.round, ] ) > 0 )
      {
        vec.phi.psi.postfusion.trans.proline <- c( vec.phi.psi.postfusion.trans.proline,
                                                   1 )
      } else
      {
        vec.phi.psi.postfusion.trans.proline <- c( vec.phi.psi.postfusion.trans.proline,
                                                   0 )
      }
    }
  }
  
  # Add column to data frame
  df.chain.resno.resid.postfusion.phi.psi.angle <- cbind( df.chain.resno.resid.postfusion.phi.psi,
                                                          vec.phi.psi.postfusion.trans.proline )
  colnames( df.chain.resno.resid.postfusion.phi.psi.angle ) <- c( colnames( df.chain.resno.resid.postfusion.phi.psi ),
                                                                  "angle_trans_proline_postfusion" )
  
  #==================================================
  # Secondary structure
  
  # Add extra column to add helix conformations
  df.chain.resno.resid.postfusion.phi.psi.angle.helix <- cbind( df.chain.resno.resid.postfusion.phi.psi.angle,
                                                                rep( 0,
                                                                     nrow( df.chain.resno.resid.postfusion.phi.psi.angle ) ),
                                                                rep( 0,
                                                                     nrow( df.chain.resno.resid.postfusion.phi.psi.angle ) ),
                                                                rep( 0,
                                                                     nrow( df.chain.resno.resid.postfusion.phi.psi.angle ) ) )
  colnames( df.chain.resno.resid.postfusion.phi.psi.angle.helix ) <- c( colnames( df.chain.resno.resid.postfusion.phi.psi.angle ),
                                                                        "helix_postfusion",
                                                                        "helix_termini_postfusion",
                                                                        "DSSP_postfusion" )
  
  
  vec.ind.subset.chain <- as.numeric( which( list.ss.dssp.postfusion$helix$chain == postfusion.chain ) )
  
  vec.helix.start <- list.ss.dssp.postfusion$helix$start[ vec.ind.subset.chain ]
  vec.helix.end <- list.ss.dssp.postfusion$helix$end[ vec.ind.subset.chain ]
  vec.index.helix <- NULL
  if( length( vec.helix.start ) > 0 )
  {
    for( i in 1:length( vec.helix.start ) )
    {
      vec.index.helix <- c( vec.index.helix,
                            vec.helix.start[ i ]:vec.helix.end[ i ] )
    }
    
    vec.resno.helix <- df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% as.character( vec.index.helix ), ]$resno_postfusion
    # Add to data frame
    df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% vec.resno.helix, ]$helix_postfusion <- 1
    # df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% c( vec.helix.start, vec.helix.end ), ]$helix_termini_postfusion <- 1
    df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$resno_postfusion %in% as.character( sort( as.vector( c( vec.helix.start, vec.helix.end ) ) ) ), ]$helix_termini_postfusion <- 1
  }
  
  
  # Add DSSP secondary structure to data frame
  vec.dssp <- list.ss.dssp.postfusion$sse
  vec.dssp[ vec.dssp == " " ] = "L"
  df.chain.resno.resid.postfusion.phi.psi.angle.helix$DSSP_postfusion <- vec.dssp
  
  
  #==================================================
  # Run pymol to mutate all positions
  print( "Run pymol to mutate all positions in the postfusion conformation" )
  vec.postfusion.no.clashes <- NULL
  df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain <- df.chain.resno.resid.postfusion.phi.psi.angle.helix[ df.chain.resno.resid.postfusion.phi.psi.angle.helix$chain_postfusion == postfusion.chain, ]
  # for( i in 1:nrow( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ) )
  vec.postfusion.no.clashes <- foreach( i = 1:nrow( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ), .combine = c ) %dopar%
  {
    # start.time <- proc.time()
    # print( paste( i,
    #               nrow( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ),
    #               sep = "/" ) )
    
    chain <- df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain[ i, 1 ]
    resno <- df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain[ i, 2 ]
    selection <- paste( chain,
                        "/",
                        resno,
                        "/",
                        sep = "" )
    
    file.pdb.postfusion.mutation <- paste( unlist( strsplit( file.pdb.postfusion,
                                                             "\\.pdb" ) ),
                                           "_",
                                           chain,
                                           resno,
                                           ".pdb",
                                           sep = "" )
    
    #==================================================
    # Generate mutated PDB file
    command <- paste( "pymol",
                      "-qc",
                      # "~/Documents/work/vrc/sbis_projects/proline_predictor/scripts/mutate.py",
                      # "mutate.py",
                      file.exe.mutate,
                      "--",
                      file.pdb.postfusion,
                      file.pdb.postfusion.mutation,
                      selection,
                      "PRO",
                      distance_angstroem )
    # Run python (PyMol) command
    pymol.out <- system( command,
                         intern = TRUE )
    if( sum( grepl( "Error",
                    pymol.out ) ) > 0 )
    {
      # vec.postfusion.no.clashes <- c( vec.postfusion.no.clashes,
      #                                 0 )
      # next
      
      return( 0 )
    }
    
    #==================================================
    # Analyze mutated protein structure
    
    # Load mutated structure
    pdb.mutant <- read.pdb( file.pdb.postfusion.mutation )
    pdb.mutant.atom <- pdb.mutant$atom
    
    # Select mutated residue
    # Do we have an insertion?
    if( is.letter( resno ) )
    {
      resno.new <- NULL
      insert.new <- NULL
      for( char in unlist( strsplit( resno,
                                     "" ) ) )
      {
        if( is.letter( char ) )
        {
          insert.new <- c( insert.new,
                           char )
        } else
        {
          resno.new <- c( resno.new,
                          char )
        }
      }
      
      resno.new <- as.numeric( paste( resno.new,
                                      collapse = "" ) )
      insert.new <- paste( insert.new,
                           collapse = "" )
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno.new &
                                                     pdb.mutant.atom$insert == insert.new &
                                                     !( is.na( pdb.mutant.atom$insert ) ), ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno.new &
                                                            pdb.mutant.atom$insert == insert.new &
                                                            !( is.na( pdb.mutant.atom$insert ) ) ), ]
      
    } else
    {
      pdb.mutant.atom.mutation <- pdb.mutant.atom[ pdb.mutant.atom$chain == chain &
                                                     pdb.mutant.atom$resno == resno, ]
      pdb.mutant.atom.non.mutation <- pdb.mutant.atom[ !( pdb.mutant.atom$chain == chain &
                                                            pdb.mutant.atom$resno == resno ), ]
    }
    # Select xyz coordinates for non-backbone atoms of mutated residue 
    pdb.mutant.atom.mutation.no.backbone <- pdb.mutant.atom.mutation[ which( !( pdb.mutant.atom.mutation$elety %in% c( "N",
                                                                                                                       "CA",
                                                                                                                       "C",
                                                                                                                       "O" ) ) ), ]
    pdb.mutant.atom.mutation.no.backbone.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.mutation.no.backbone[ , 9:11 ] ) ) )
    
    # Select xyz coordinates for all other residues
    pdb.mutant.atom.non.mutation.xyz <- as.vector( t( as.matrix( pdb.mutant.atom.non.mutation[ , 9:11 ] ) ) )
    
    # Calculate distance matrix
    mat.dist <- dist.xyz( pdb.mutant.atom.mutation.no.backbone.xyz,
                          pdb.mutant.atom.non.mutation.xyz )
    
    # Generate vectors representing atom types for mutated residue and rest
    vec.atom.type.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.mutation.no.backbone$elety ) )
    {
      elety <- pdb.mutant.atom.mutation.no.backbone$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.mutation <- c( vec.atom.type.mutation,
                                   elety.split[ which( unlist( lapply( elety.split,
                                                                       is.letter ) ) )[ 1 ] ] )
    }
    
    vec.atom.type.non.mutation <- NULL
    for( j in 1:length( pdb.mutant.atom.non.mutation$elety ) )
    {
      elety <- pdb.mutant.atom.non.mutation$elety[ j ]
      elety.split <- unlist( strsplit( elety,
                                       "" ) )
      vec.atom.type.non.mutation <- c( vec.atom.type.non.mutation,
                                       elety.split[ which( unlist( lapply( elety.split,
                                                                           is.letter ) ) )[ 1 ] ] )
    }
    
    # Calculate MolProbity clash constraint
    df.atom.type.mutation.nonmutation <- data.frame( expand.grid( unique( vec.atom.type.mutation ),
                                                                  unique( vec.atom.type.non.mutation ) ) )
    df.atom.type.mutation.nonmutation[ , 1 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 1 ] )
    df.atom.type.mutation.nonmutation[ , 2 ] <- as.vector( df.atom.type.mutation.nonmutation[ , 2 ] )
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- NULL
    for( j in 1:nrow( df.atom.type.mutation.nonmutation ) )
    {
      r1 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 1 ], ]$rvdw
      r2 <- elements[ elements$symb == df.atom.type.mutation.nonmutation[ j, 2 ], ]$rvdw
      df.atom.type.mutation.nonmutation.r1.r2.0.9 <- rbind( df.atom.type.mutation.nonmutation.r1.r2.0.9,
                                                            c( df.atom.type.mutation.nonmutation[ j, ],
                                                               r1 + r2 - overlap.cutoff ) )
    }
    df.atom.type.mutation.nonmutation.r1.r2.0.9 <- data.frame( df.atom.type.mutation.nonmutation.r1.r2.0.9 )
    colnames( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) <- c( "mutation",
                                                                  "nonmutation",
                                                                  "value" )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$mutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation <- as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$nonmutation ) )
    df.atom.type.mutation.nonmutation.r1.r2.0.9$value <- as.numeric( as.vector( unlist( df.atom.type.mutation.nonmutation.r1.r2.0.9$value ) ) )
    
    mat.r1.r2.0.9 <- matrix( 0,
                             length( vec.atom.type.mutation ),
                             length( vec.atom.type.non.mutation ) )
    
    for( j in 1:nrow( df.atom.type.mutation.nonmutation.r1.r2.0.9 ) )
    {
      vec.i <- which( vec.atom.type.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$mutation )
      vec.j <- which( vec.atom.type.non.mutation == df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$nonmutation )
      
      mat.r1.r2.0.9[ vec.i, vec.j ] <- df.atom.type.mutation.nonmutation.r1.r2.0.9[ j, ]$value
    }
    
    # Get number of clashes
    mat.clash <- mat.dist < mat.r1.r2.0.9
    
    # vec.postfusion.no.clashes <- c( vec.postfusion.no.clashes,
    #                                 sum( mat.clash ) )
    # print( vec.postfusion.no.clashes )
    # print( ( proc.time() - start.time )[ 3 ] )
    sum( mat.clash )
  }
  # vec.postfusion.no.clashes <- unlist( vec.postfusion.no.clashes )
  
  df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes <- cbind( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain,
                                                                                vec.postfusion.no.clashes )
  colnames( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes ) <- c( colnames( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain ),
                                                                                        "no_clash_postfusion" )
  write.table( df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes,
               file = paste( output_prefix,
                             "_postfusion.txt",
                             sep = "" ) )
}





#==================================================
#==================================================
#==================================================
# Generate output files
if( option == 1 )
{
  #==================================================
  # Merge two data frames
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes <- merge( df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes,
                                                                                    df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes,
                                                                                    by.x = "resno_prefusion",
                                                                                    by.y = "resno_postfusion",
                                                                                    all = TRUE )
  
  vec.dist.phi.prefusion.postfusion <- abs( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes$phi_prefusion - df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes$phi_postfusion )
  vec.dist.phi.prefusion.postfusion[ which( vec.dist.phi.prefusion.postfusion > 180 ) ] <- 360 - vec.dist.phi.prefusion.postfusion[ which( vec.dist.phi.prefusion.postfusion > 180 ) ]
  vec.dist.psi.prefusion.postfusion <- abs( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes$psi_prefusion - df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes$psi_postfusion )
  vec.dist.psi.prefusion.postfusion[ which( vec.dist.psi.prefusion.postfusion > 180 ) ] <- 360 - vec.dist.psi.prefusion.postfusion[ which( vec.dist.psi.prefusion.postfusion > 180 ) ]
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles <- cbind( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes,
                                                                                               vec.dist.phi.prefusion.postfusion,
                                                                                               vec.dist.psi.prefusion.postfusion )
  colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ) <- c( colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes ),
                                                                                                       "dist_phi_prefusion_postfusion",
                                                                                                       "dist_psi_prefusion_postfusion" )
  
  #==================================================
  # Could you add an additional column called “compatible with prefusion”,  1 if 1,0,0 is at column F-H, 0 otherwise?
  vec.compatible.with.prefusion <- NULL
  for( i in 1:nrow( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ) )
  {
    if( is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_prefusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_prefusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_prefusion[ i ] ) )
    {
      vec.compatible.with.prefusion <- c( vec.compatible.with.prefusion,
                                          NA )
    } else if( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_prefusion[ i ] == 1 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_prefusion[ i ] == 0 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_prefusion[ i ] == 0 )
    {
      vec.compatible.with.prefusion <- c( vec.compatible.with.prefusion,
                                          1 )
    } else
    {
      vec.compatible.with.prefusion <- c( vec.compatible.with.prefusion,
                                          0 )
    }
  }
  
  #==================================================
  # Could you add an additional column called “compatible with postfusion”,  0 if 1,0,0 is at column M-O, 1 otherwise?
  vec.compatible.with.postfusion <- NULL
  for( i in 1:nrow( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ) )
  {
    if( is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_postfusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_postfusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_postfusion[ i ] ) )
    {
      vec.compatible.with.postfusion <- c( vec.compatible.with.postfusion,
                                           NA )
    } else if( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_postfusion[ i ] == 1 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_postfusion[ i ] == 0 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_postfusion[ i ] == 0 )
    {
      vec.compatible.with.postfusion <- c( vec.compatible.with.postfusion,
                                           1 )
    } else
    {
      vec.compatible.with.postfusion <- c( vec.compatible.with.postfusion,
                                           0 )
    }
  }
  
  #==================================================
  # Compatible prefusion and not postfusion
  vec.compatible.with.prefusion.not.postfusion <- NULL
  for( i in 1:length( vec.compatible.with.prefusion ) )
  {
    if( is.na( vec.compatible.with.prefusion[ i ] ) |
        is.na( vec.compatible.with.postfusion[ i ] ) )
    {
      vec.compatible.with.prefusion.not.postfusion <- c( vec.compatible.with.prefusion.not.postfusion,
                                                         NA )
    } else if( vec.compatible.with.prefusion[ i ] == 1 &
               vec.compatible.with.postfusion[ i ] == 0 )
    {
      vec.compatible.with.prefusion.not.postfusion <- c( vec.compatible.with.prefusion.not.postfusion,
                                                         1 )
    } else
    {
      vec.compatible.with.prefusion.not.postfusion <- c( vec.compatible.with.prefusion.not.postfusion,
                                                         0 )
    }
  }
  
  
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC <- cbind( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles,
                                                                                                           vec.compatible.with.prefusion,
                                                                                                           vec.compatible.with.postfusion,
                                                                                                           vec.compatible.with.prefusion.not.postfusion )
  colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC ) <- c( colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ),
                                                                                                                   "compatiblePrefusion",
                                                                                                                   "compatiblePostfusion",
                                                                                                                   "compatiblePrefusionNotPostfusion" )
  
  
  #==================================================
  #==================================================
  # Save data frame
  write.csv( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC,
             file = paste( output_prefix,
                           "_DF.csv",
                           sep = "" ),
             row.names = FALSE )
} else if( option == 2 )
{
  #==================================================
  # prepare data frame
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles <- df.chain.resno.resid.prefusion.phi.psi.angle.helix.chain.noclashes
  
  
  #==================================================
  # Could you add an additional column called “compatible with prefusion”,  1 if 1,0,0 is at column F-H, 0 otherwise?
  vec.compatible.with.prefusion <- NULL
  for( i in 1:nrow( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ) )
  {
    if( is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_prefusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_prefusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_prefusion[ i ] ) )
    {
      vec.compatible.with.prefusion <- c( vec.compatible.with.prefusion,
                                          NA )
    } else if( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_prefusion[ i ] == 1 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_prefusion[ i ] == 0 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_prefusion[ i ] == 0 )
    {
      vec.compatible.with.prefusion <- c( vec.compatible.with.prefusion,
                                          1 )
    } else
    {
      vec.compatible.with.prefusion <- c( vec.compatible.with.prefusion,
                                          0 )
    }
  }
  
  #==================================================
  # Add compatibility information
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC <- cbind( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles,
                                                                                                           vec.compatible.with.prefusion )
  colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC ) <- c( colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ),
                                                                                                                   "compatiblePrefusion" )
  
  #==================================================
  #==================================================
  # Save data frame
  write.csv( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC,
             file = paste( output_prefix,
                           "_DF.csv",
                           sep = "" ),
             row.names = FALSE )
} else if( option == 3 )
{
  #==================================================
  # prepare data frame
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles <- df.chain.resno.resid.postfusion.phi.psi.angle.helix.chain.noclashes
  
  #==================================================
  # Could you add an additional column called “compatible with postfusion”,  0 if 1,0,0 is at column M-O, 1 otherwise?
  vec.compatible.with.postfusion <- NULL
  for( i in 1:nrow( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ) )
  {
    if( is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_postfusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_postfusion[ i ] ) |
        is.na( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_postfusion[ i ] ) )
    {
      vec.compatible.with.postfusion <- c( vec.compatible.with.postfusion,
                                           NA )
    } else if( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$angle_trans_proline_postfusion[ i ] == 1 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$helix_postfusion[ i ] == 0 &
               df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles$no_clash_postfusion[ i ] == 0 )
    {
      vec.compatible.with.postfusion <- c( vec.compatible.with.postfusion,
                                           1 )
    } else
    {
      vec.compatible.with.postfusion <- c( vec.compatible.with.postfusion,
                                           0 )
    }
  }
  
  #==================================================
  # Add compatibility information
  df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC <- cbind( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles,
                                                                                                           vec.compatible.with.postfusion )
  colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC ) <- c( colnames( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles ),
                                                                                                                   "compatiblePostfusion" )
  
  #==================================================
  #==================================================
  # Save data frame
  write.csv( df.chain.resno.resid.prefusion.postfusion.phi.psi.angle.helix.noclashes.distangles.criteriaGYC,
             file = paste( output_prefix,
                           "_DF.csv",
                           sep = "" ),
             row.names = FALSE )
}





#==================================================
# Measure and print out time needed for the script
end.main <- proc.time()
duration.main <- end.main-start.main
print( paste( "Script duration:", round( duration.main[3] / 60, 2 ), "min") )
#==================================================

# Main end
#==================================================
#==================================================
#==================================================