setwd("/Users/xji3/GitFolders/EDN_ECP/SimulationStudy/")
# Now read in HMM results from plots
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/Tract_", toString(tract), '.0_HKY/sim_', 
                      toString(sim), '/HMM_EDN_ECP_HKY_rv_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}

# Now read in HMM results from plots for HalfTau case
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/HalfTau/Tract_", toString(tract), '.0_HKY/sim_', 
                      toString(sim), '/HMM_EDN_ECP_HKY_rv_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("HalfTau_HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}

# Now read in HMM results from plots for TenthTau case
Tract.list <- c(3.0, 10.0, 50.0, 100.0, 200.0, 300.0, 400.0, 500.0)
for(tract in Tract.list){
  hmm.tract.plots <- NULL
  for(sim in 1:100){
    hmm.plot <- paste("./plot/TenthTau/Tract_", toString(tract), '.0_HKY/sim_', 
                      toString(sim), '/HMM_EDN_ECP_HKY_rv_lnL_sim_',
                      toString(sim), '_1D_surface.txt', sep = "")
    if (file.exists(hmm.plot)){
      lnL.surface <- read.table(hmm.plot)
      max.idx <- which.max(lnL.surface[, 2])
      new.summary <- matrix(c(3.0*exp(-lnL.surface[max.idx, 1]), lnL.surface[max.idx, 2]), 2, 1)
      rownames(new.summary) <- c("tract in nt", "lnL")
      colnames(new.summary) <- paste("sim_", toString(sim), sep = "")
      hmm.tract.plots <- cbind(hmm.tract.plots, new.summary)     
    }
  }
  assign(paste("TenthTau_HMM_Tract_", toString(tract), "_plot", sep = ""), hmm.tract.plots)
}