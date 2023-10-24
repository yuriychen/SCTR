#translation related functions
condon_t = list(
  'TTT'='F', 'TTC'='F', 'TTA'='L', 'TTG'='L',
  'TCT'='S', 'TCC'='S', 'TCA'='S', 'TCG'='S',
  'TAT'='Y', 'TAC'='Y', 'TAA'='*', 'TAG'='*',
  'TGT'='C', 'TGC'='C', 'TGA'='*', 'TGG'='W',
  'CTT'='L', 'CTC'='L', 'CTA'='L', 'CTG'='L',
  'CCT'='P', 'CCC'='P', 'CCA'='P', 'CCG'='P',
  'CAT'='H', 'CAC'='H', 'CAA'='Q', 'CAG'='Q',
  'CGT'='R', 'CGC'='R', 'CGA'='R', 'CGG'='R',
  'ATT'='I', 'ATC'='I', 'ATA'='I', 'ATG'='M',
  'ACT'='T', 'ACC'='T', 'ACA'='T', 'ACG'='T',
  'AAT'='N', 'AAC'='N', 'AAA'='K', 'AAG'='K',
  'AGT'='S', 'AGC'='S', 'AGA'='R', 'AGG'='R',
  'GTT'='V', 'GTC'='V', 'GTA'='V', 'GTG'='V',
  'GCT'='A', 'GCC'='A', 'GCA'='A', 'GCG'='A',
  'GAT'='D', 'GAC'='D', 'GAA'='E', 'GAG'='E',
  'GGT'='G', 'GGC'='G', 'GGA'='G', 'GGG'='G'
)

reverse_codon_t = list(
  'A'= c('GCT','GCC','GCA','GCG'),
  'C'= c('TGT','TGC'),
  'D'= c('GAT','GAC'),
  'E'= c('GAA','GAG'),
  'F'= c('TTT','TTC'),
  'G'= c('GGT','GGC','GGA','GGG'),
  'H'= c('CAT','CAC'),
  'I'= c('ATT','ATC','ATA'),
  'K'= c('AAA','AAG'),
  'L'= c('TTA','TTG','CTT','CTC','CTA','CTG'),
  'M'= c('ATG'),
  'N'= c('AAT','AAC'),
  'P'= c('CCT','CCC','CCA','CCG'),
  'Q'= c('CAA','CAG'),
  'R'= c('CGT','CGC','CGA','CGG','AGA','AGG'),
  'S'= c('TCT','TCC','TCA','AGT','AGC'),
  'T'= c('ACT','ACC','ACA','ACG'),
  'V'= c('GTT','GTC','GTA','GTG'),
  'W'= c('TGG'),
  'Y'= c('TAT','TAC'),
  '*'= c('TAG','TGA','TAA'),
  'U'= c('TCG')
)

nc_to_aa = function(codon_nc, codon_table=condon_t){
  codon_nc <- toupper(codon_nc)
  return(as.character(codon_table[codon_nc]))
}

ncseq_to_aaseq = function(nc_seq, codon_table=codon_t){
  aa_seq = ''
  for (i in seq(1,nchar(nc_seq),3)) {
    aa_seq = paste(aa_seq, as.character(nc_to_aa(substr(nc_seq,i,i+2))), sep='')
  }
  return(aa_seq)
}

aa_to_nc = function(aa, reverse_codon_table=reverse_codon_t){
  aa <- toupper(aa)
  return(reverse_codon_table[[aa]])
}

#weights change here
###
weights = read.table('weights_ml_20230613.txt',sep=' ',header=FALSE)
colnames(weights) = c('A','T','C','G')

#score calculation related functions
matrix_construct_adjust = function(seq, contain_mid=TRUE, oneside_num=6){
  seq = toupper(seq)
  
  if (contain_mid) {
    if (nchar(seq) != (oneside_num *2 + 3)){
      stop('seq len is wrong!')
    }
  }
  else{
    if (nchar(seq) != (oneside_num * 2)){
      stop('seq len is wrong!')
    }
  }
  
  if (contain_mid) {
    seq = paste(substr(seq,1,oneside_num), substr(seq, oneside_num + 4,oneside_num * 2 + 3), sep='')
  }
  
  seq_matrix = data.frame(matrix(0, nrow=nchar(seq), ncol=4))
  colnames(seq_matrix) <- c('A','T','C','G')
  
  for (i in seq(1,nchar(seq))){
    c = substr(seq,i,i)
    seq_matrix[i,c] = 1
  }
  
  return(seq_matrix)
}

score_cal_adjust = function(seq_matrix, weights_matrix=weights, oneside_num=6){
  weights_matrix = weights_matrix[((30 - oneside_num + 1) : (oneside_num + 30)) ,]
  result_matrix = weights_matrix * seq_matrix
  result = (sum(result_matrix) - sum(apply(weights_matrix, 1, min)))/(sum(apply(weights_matrix, 1, max)) - sum(apply(weights_matrix, 1, min)))
  return(result)
}

seq_to_score_matrix = function(seq, seq_cutoff=6, oneside_num=6, weights_matrix=weights, contain_mid=TRUE){
  if (contain_mid) {
    if (nchar(seq) < (oneside_num *2 + 3)){
      stop('seq len is too short!')
    }
  }
  
  nc_site_list = c()
  aa_site_list = c()
  site_translate_list = c()
  subseq_list = c()
  score_list = c()
  
  seq = toupper(seq)
  
  for (index in seq(1,(nchar(seq)/3))){
    if ((index * 3 - 2) < seq_cutoff){
      next
    }
    if ((index * 3) > (nchar(seq) - seq_cutoff)){
      break
    }
    
    nc_site_list = append(nc_site_list, (index * 3 - 2))
    aa_site_list = append(aa_site_list, (index))
    site_translate_list = append(site_translate_list, as.character(nc_to_aa(substr(seq,(index * 3 - 2),(index * 3)))))
    
    subseq = substr(seq, (index * 3 - 2 - oneside_num), (index * 3 + oneside_num))
    subseq_list = append(subseq_list, subseq)
    score_list = append(score_list, score_cal_adjust(matrix_construct_adjust(subseq, contain_mid, oneside_num), weights_matrix, oneside_num))
  }
  
  result_df = data.frame(nc_site_list, aa_site_list, site_translate_list, subseq_list, score_list)
  colnames(result_df) = c('NC_site','AA_site','Site_aa','Subseq','Score')
  
  return(result_df)
}

#nc sequence redesign functions
sequence_list_generate = function(aa_seq,reverse_codon_table=reverse_codon_t,contain_mid=TRUE,oneside_num=6){
  aa_seq = toupper(aa_seq)
  
  if (contain_mid) {
    if (nchar(aa_seq) != ((oneside_num *2 + 3) / 3)){
      stop('seq len is wrong!')
    }
  }
  else{
    if (nchar(aa_seq) != ((oneside_num * 2) / 3)){
      stop('seq len is wrong!')
    }
  }
  
  if (contain_mid == FALSE) {
    aa_seq = paste(substr(aa_seq,1,(oneside_num / 3)),'U', substr(aa_seq, oneside_num / 3 + 1, (oneside_num * 2) / 3), sep='')
  }
  
  seq_list = c()
  aa_list = sapply(seq(from=1, to=nchar(aa_seq), by=1), function(i) substr(aa_seq, i, i))
  for (aa in aa_list) {
    codon_list = aa_to_nc(aa)
    
    if (length(seq_list) == 0) {
      seq_list = codon_list
    }
    else{
      temp_list = c()
      for (upseq in seq_list){
        for (codon in codon_list){
          temp_list = append(temp_list, paste(upseq,codon,sep=''))
        }
      }
      seq_list = temp_list
    }
  }
  return(seq_list)
}

sequence_list_score_calc = function(seq_list,weights_matrix=weights,contain_mid=TRUE,oneside_num=6){
  score_list = c()
  for (s in seq_list){
    score = score_cal_adjust(matrix_construct_adjust(s,contain_mid,oneside_num),weights_matrix,oneside_num)
    score_list = append(score_list,score)
  }
  return(score_list)
}

aa_seq_score_matrix = function(aa_seq,weights_matrix=weights,reverse_codon_table=reverse_codon_t,contain_mid=TRUE,oneside_num=6){
  seqs = sequence_list_generate(aa_seq,reverse_codon_table,contain_mid)
  scores = sequence_list_score_calc(seqs,weights_matrix,contain_mid,oneside_num)
  score_matrix = data.frame(seqs,scores)
  colnames(score_matrix) = c('Seq','Score')
  score_matrix = score_matrix[order(score_matrix$Score,decreasing = TRUE),]
  rownames(score_matrix) = seq(1,nrow(score_matrix))
  score_matrix$Seq = factor(score_matrix$Seq, levels=score_matrix$Seq)
  return(score_matrix)
}

ncseq_redesign_score_matrix = function(ncseq, weights_matrix=weights,contain_mid=TRUE,oneside_num=6){
  if (nchar(ncseq) != ((oneside_num *2 + 3))){
    stop('seq len is wrong!')
  }
  
  ncseq = toupper(ncseq)
  aa_seq = ncseq_to_aaseq(ncseq)[[1]]
  aa_ori = aa_seq
  aa_seq = paste(substr(aa_seq,1,(oneside_num / 3)),'U',substr(aa_seq, oneside_num / 3 + 2, (oneside_num * 2 + 3) / 3), sep='')
  score_ori = score_cal_adjust(matrix_construct_adjust(ncseq,contain_mid,oneside_num),weights_matrix,oneside_num)
  score_df = aa_seq_score_matrix(aa_seq,weights_matrix,reverse_codon_t,contain_mid,oneside_num)  
  score_df$delta_score = score_df$Score - score_ori
  score_df$Score_original = score_ori
  score_df$AA = aa_seq
  score_df$AA_original = aa_ori
  
  return(score_df)
}

#aa design functions
aaseq_redesign_score_matrix = function(aaseq, weights_matrix=weights,contain_mid=TRUE,oneside_num=6){
  if (nchar(aaseq) != ((oneside_num * 2 + 3) / 3)){
    stop('seq len is wrong!')
  }
  
  aaseq = toupper(aaseq)
  aa_ori = aaseq
  aaseq = paste(substr(aaseq,1,(oneside_num / 3)),'U',substr(aaseq, oneside_num / 3 + 2, (oneside_num * 2 + 3) / 3), sep='')
  score_df = aa_seq_score_matrix(aaseq,weights_matrix,reverse_codon_t,contain_mid,oneside_num)  
  score_df$AA = aaseq
  score_df$AA_original = aa_ori
  
  return(score_df)
}







