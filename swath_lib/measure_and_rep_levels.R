#
ms.rep <- c('run_id', 'run_code', 'ms')
tech.rep <- c(ms.rep, 'sample_code', 'tech_id', 'tech')
bio.rep <- c(tech.rep, 'bio', 'state')
frg.lev <- c('PeakName', 'fragment_id', 'fragment_series', 'fragment_number', 'fragment_charge', 'clean', 'cluster')
prec.lev <- c(frg.lev, 'precursor_id', 'precursor_charge')
pep.lev <- c(prec.lev, 'peptide_id', 'clean_peptide_id', 'peptide_sequence', 'clean_peptide_sequence', 'modifications', 'uniq', 'tryptic', 'selected_modifications')
