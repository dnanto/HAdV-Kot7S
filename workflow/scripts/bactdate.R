#!/usr/bin/env Rscript --vanilla

main <- function(snakemake) {
  library(tidyverse)
  library(BactDating)
  
  path.in.tree <- snakemake@input[["tree"]]
  path.in.node <- snakemake@input[["node"]]
  path.in.stat <- snakemake@input[["stat"]]
  path.in.meta <- snakemake@input[["meta"]]
  path.in.fast <- snakemake@input[["fast"]]
  path.out.tree <- snakemake@output[["tree"]]
  path.out.json <- snakemake@output[["json"]]
  path.out.tsv <- snakemake@output[["tsv"]]
  params <- snakemake@params
  threads <- snakemake@threads

  stopifnot(file.exists(path.in.tree, path.in.node, path.in.stat, path.in.meta))

  # seeds
  set.seed(params$seed)
  seeds <- runif(max = .Machine$integer.max, n = length(params$models) * params$reps)

  # simulate
  df.meta <- read_tsv(path.in.meta, col_types = cols(.default = "c"))
  phy <- loadGubbins(stringr::str_remove(path.in.tree, ".final_tree.tre$"))
  phy$tip.date <- (
    with(df.meta, date[match(phy$tip.label, accession)]) %>%
      lubridate::parse_date_time(c("%Y-%m-%d", "%Y-%m", "%Y")) %>%
      lubridate::decimal_date()
  )
  phy <- initRoot(phy, phy$tip.date, useRec = T)
  root <- dirname(path.out.tree)
  paths <- with(
    expand.grid(model = params$models, rep = 1:params$reps),
    parallel::mcmapply(
      function(model, rep, seed) {
        path <- file.path(root, paste0("bac-", model, "-", rep, ".qs"))
        if (!file.exists(path)) {
          set.seed(seed)
          run <- bactdate(phy, phy$tip.date, nbIts = params$nbIts, thin = params$thin, model = model, useRec = T)
          run$seed <- seed
          run$rep <- rep
          # output run
          qs::qsave(run, path)
        }
        path
      },
      model = model, rep = rep, seed = seeds, mc.cores = threads
    )
  )

  # model comparison
  key <- c("likelihood", "mu", "sigma", "alpha", "prior")
  aln <- ncol(ape::read.dna(path.in.fast, format = "fasta", as.character = T))
  df.diagnostic <- (
    parallel::mclapply(
      paths,
      function(path) {
        run <- qs::qread(path)
        key <- c("likelihood", "mu", "sigma", "alpha", "prior")
        rec <- with(run, record[max(1, round(nrow(record) * params$burn)):nrow(record), key])
        est <- (
          apply(rec[, key], 2, summary) %>%
            t() %>%
            apply(1, function(row) sprintf("%.2e [%.2e;%.2e]", row["Median"], row["1st Qu."], row["3rd Qu."])) %>%
            setNames(paste0("est.", key))
        )
        ess <- setNames(coda::effectiveSize(rec), paste0("ess.", key))
        rate <- ((
          enframe(capture.output(run)[14], NULL) %>%
            separate(value, c("_", "rate", "__", "low", "high", "___"), "[=\\[\\]; ]") %>%
            with(c(as.numeric(rate) / aln, as.numeric(low) / aln, as.numeric(high) / aln)) %>%
            formatC(format = "e", digits = 2) %>%
            {
              paste(.[1], " [", .[2], ";", .[3], "]", sep = "")
            })
        )
        c(
          path = path,
          model = as.character(run$model),
          rep = run$rep,
          seed = run$seed,
          dic = run$dic,
          rootprob = run$rootprob,
          est.rate = rate,
          est.root_date_likeliest = str_split_fixed(capture.output(run)[18], "=", 2)[, 2],
          est,
          ess
        )
      },
      mc.cores = threads
    ) %>%
      bind_rows() %>%
      mutate(
        across(starts_with("ess."), as.numeric),
        ess.threshold = params$ess,
        pass = (
          (params$ess <= ess.sigma | model == "poisson" | model == "strictgamma") +
            (params$ess <= ess.likelihood) +
            (params$ess <= ess.mu) +
            (params$ess <= ess.alpha) +
            (params$ess <= ess.prior)
        )
      ) %>%
      arrange(desc(pass), desc(dic), across(starts_with("ess."), desc))
  )

  # output
  # reload best model
  run <- (
    pull(df.diagnostic, path) %>%
      head(1) %>%
      qs::qread(nthreads = threads)
  )
  # output tree
  ape::write.tree(run$tree, path.out.tree)
  # output json
  root <- min(run$tree$edge[, 1])
  iden <- c(root, run$tree$edge[, 2])
  label <- c(run$tree$tip.label, run$tree$node.label)[iden]
  data <- data.frame(
    branch_length = c(0, run$tree$edge.length),
    clock_length = c(0, run$tree$edge.length),
    mutation_length = c(0, run$tree$subs),
    num_date = c(leafDates(run$tree), nodeDates(run$tree))[iden],
    num_date_confidence = I(run$CI[iden, ])
  )
  write_lines(rjson::toJSON(list(nodes = split(data, label)), indent = 1), path.out.json)
  # output tsv
  write_tsv(select(df.diagnostic, -path), path.out.tsv)
}

capture.output(
  capture.output(main(snakemake), file = snakemake@log[["out"]], type = "output"),
  file = snakemake@log[["msg"]], type = "message"
)
