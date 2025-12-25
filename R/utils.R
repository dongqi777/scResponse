#-----------------------------------------------------------------------------------------------------
#这个函数不太重要，训练时使用的
quintile_normalize <- function(mat) {
  quintiles <- apply(mat, 2, function(x) {  quantile(x, probs = seq(0, 1, 0.2), na.rm = TRUE)})   #  计算每个样本的五分位数
  ref_dist <- rowMeans(quintiles, na.rm = TRUE)                                                   #  计算参考分布（所有样本的平均五分位数）

  ##########  创建归一化矩阵
  normalized <- apply(mat, 2, function(col)  {  # col=mat[,2]  #test
    q <- quantile(col, probs = seq(0, 1, 0.2), na.rm = TRUE)                    #  计算当前样本的五分位数
    approx_fun <- approxfun(q, ref_dist, rule = 2, method = "linear")
    #approxfun  线性插值函数  ,根据给定的数据点 (x,y) 生成一个近似函数，可以用来计算任意新x 值对应的y 值。它基于线性插值（默认）或样条插值，适用于数据填充、平滑和预测。  # rule:外推规则：1(默认，返回NA)或2(使用端点值)。
    approx_fun(col)  }  )  # 应用插值
  dimnames(normalized) <- dimnames(mat)  # 保持原始行名和列名
  return(normalized)   }

#-----------------------------------------------------------------------------------------------------
#添加函数，解决过度依赖
build_rankings_simple <- function(mat) {
  # mat: gene x cell
  apply(mat, 2, function(x) {
    rank(-x, ties.method = "average")
  })
}

#-----------------------------------------------------------------------------------------------------

calculate_deviation <- function(gene_Diff, nperm=3000, gene_list, topN ) {  ### gene_list:表达谱所有基因
  # gene_list should be decreasingly ordered based on their expression
  # assuming the top genes in each sample are usually within the top half highly expressed genes across samples.
  # topN=topN
  if(topN > round(length(gene_list)/2) ) {
    gene_list <- gene_list[1: (2*round(length(gene_list)/3) )  ]  } else {
      gene_list <- gene_list[1:round(length(gene_list)/2)] }  #选取高表达前1/2的基因。 gene_list:所有基因表达均值。round,四舍五入

  #GSEA的目的是判断gene_Diff里面的成员在 list L(topN表达基因)里面是随机分布还是主要聚集在L的顶部或底部。
  ES_null <- numeric(nperm)  ## 生成nperm个0值
  for (i in 1:nperm) {  # i=1
    #1.  按顺序随机选取基因， 随机nperm次
    sorted_indices <- sort(sample.int(length(gene_list),topN))  ##在高表达前1/2的基因中随机抽取topN个基因，升序（即基因表达值降序）。 ##sample.int, 从整数范围内随机抽取样本
    topN_genes <- gene_list[sorted_indices]
    cumulative_scores <- cumsum(gene_Diff[topN_genes])          # 利用topN_genes（500个），及其Diff分数来计算累积排秩得分
    #cumsum()是R语言中用于计算累积和（cumulative sum）的函数，它会返回一个向量，其中每个元素都是输入向量从开始到当前位置所有元素的和。
    # 2. 计算最大正偏差和最大负偏差
    max_positive <- max(cumulative_scores) ;  max_negative <- min(cumulative_scores)
    # 3. 确定主导方向和其位置
    if (abs(max_positive) >= abs(max_negative)) {   #abs绝对值
      direction <- "positive"  ;  max_dev <- max_positive } else {  direction <- "negative"  ;  max_dev <- max_negative }
    # 4. 计算这次随机打乱的曲线的ES
    ES_null[i] <- max_dev }

  return(   list(random_mean = mean(ES_null),   random_sd = sd(ES_null)   )    )    }

#-----------------------------------------------------------------------------------------------------
calculate_enrichment_score <- function(cumulative_scores, random_mean=0, random_sd=1,  nGenes ) {
  #  nGenes <- nGenes   # 检测到的基因总数  #cumulative_scores每个样本/细胞对应一组cumulative_scores

  # 1. 计算最大正偏差和最大负偏差
  max_positive <- max(cumulative_scores) ;  max_negative <- min(cumulative_scores)

  # 2. 确定主导方向和其位置
  if (  abs(max_positive) >= abs(max_negative) ) {
    direction <- "positive";  max_dev <- max_positive  ;   max_pos <- which.max(cumulative_scores)   } else {
      direction <- "negative";  max_dev <- max_negative  ;   max_pos <- which.min(cumulative_scores)   }

  # 3. 计算位置权重（早期偏差更重要）
  position_weight <- (1 - (max_pos / nGenes) )^2   ##max_pos / nGenes →  ，  1 - (max_pos / nGenes) → 反转权重（使靠前的基因获得更高值），  ^2 → 平方操作，强化高权重基因的优势
  # 1 - (max_pos / nGenes)：越靠前的基因（max_pos 小），此值越接近1；越靠后的基因（max_pos 大），此值越接近0。
  # 平方操作 ^2：进一步放大靠前基因的权重，削弱靠后基因的影响（非线性衰减）。

  # 4. 计算方向强度（主导偏差与反向偏差的比例）
  direction_strength <- if ( direction == "positive" ) {
    abs(max_positive) / (abs(max_positive) + abs(max_negative)) } else {   #结果接近 1：正向信号占主导（如 max_positive 远大于 max_negative）。
      abs(max_negative) / (abs(max_positive) + abs(max_negative)) }
  #若 max_negative = 0，比值为1（完全正向信号）。若 max_positive = 0，比值为0（完全负向信号）。

  # 5. 计算标准化综合得分
  relative_dev <- (max_dev - random_mean) / random_sd                      ##标准化，即z-score。
  normalized_score <- relative_dev * position_weight * direction_strength  ## 最终
  ########################！！！！！！
  # 6. 添加符号表示方向
  #normalized_score <- if (direction == "positive") normalized_score else -normalized_score

  # 7. 计算偏差区域占比
  n <- length(cumulative_scores)     ## DQ添加
  positive_area <- sum(cumulative_scores > 0) / n   ;    negative_area <- sum(cumulative_scores < 0) / n
  return(list(  normalized_score = normalized_score, direction = direction,   max_deviation = max_dev,      max_position = max_pos,
                position_weight = position_weight,   direction_strength = direction_strength,     positive_area = positive_area,
                negative_area = negative_area,       cumulative_scores = cumulative_scores  )      )}


#-----------------------------------------------------------------------------------------------------
#这个函数是不让那些rda进入全局，rda被读入时会生成gene_giff和gene_list，用户是不知道的，如果用户全局中恰好有同名的变量，使用我们的包就会被覆盖掉，这是绝对不允许的。
.load_cancer_data_internal <- function(cancer) {

  ## 1️⃣ 从包内 data 目录获取 rda 文件名（BC / CRC / BLCA）
  data_dir <- system.file("data", package = "scResponse")
  available <- sub("\\.rda$", "", list.files(data_dir, pattern = "\\.rda$"))

  if (!cancer %in% available) {
    stop(
      sprintf(
        "Unsupported cancer type '%s'. Available types: %s",
        cancer,
        paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  ## 2️⃣ 用隔离环境加载，避免污染用户全局
  env <- new.env(parent = emptyenv())

  data(list = cancer, package = "scResponse", envir = env)

  ## 3️⃣ 校验 rda 内部结构（你的真实假设）
  if (!exists("gene_Diff", envir = env, inherits = FALSE) ||
      !exists("gene_list", envir = env, inherits = FALSE)) {
    stop(
      sprintf(
        "Dataset '%s' must contain objects: gene_Diff and gene_list",
        cancer
      ),
      call. = FALSE
    )
  }

  ## 4️⃣ 返回，不泄露到全局
  list(
    gene_Diff = env$gene_Diff,
    gene_list = env$gene_list
  )
}


