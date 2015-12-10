#!/bin/bash
# all taxa models
# switch to just the fauna models
R CMD BATCH --vanilla ../R/turnover_model.r
R CMD BATCH --vanilla ../R/turnover_process.r
