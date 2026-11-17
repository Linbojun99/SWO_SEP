# =====================================================================================================
# Shared functions for comparing the walkthrough results with the archived results and for the bootstrap samples
# 用于将本流程结果与存档结果比较以及用于自助抽样的共用函数
#
# Every step sources this file. compare_table() reads a new table and the archived table of the same name and
# prints whether they match or the largest difference. compare_figure() compares the pixels of a new PNG figure
# with the archived PNG of the same name. The other functions sort tables, check numbers and draw the bootstrap
# samples of steps 22, 25, 27 and 30.
# 每个步骤都读入本文件。compare_table() 读取一张新表及同名的存档表，输出两者是否一致或最大差值。compare_figure()
# 比较新绘 PNG 图与同名存档 PNG 图的像素。其余函数用于表格排序、数值核对，以及第 22、25、27 和 30 步的自助抽样。
#
# Reads / 读取:  the tables and figures named by each step / 各步骤指定的表和图
# Writes / 写出: nothing / 无
# =====================================================================================================

out <- "analysis_walkthrough/output"
archive <- "results/directed_fishery"

# Some archived tables store logical columns as True and False, the tables of this walkthrough as TRUE and
# FALSE; read both as logical.
# 部分存档表把逻辑列写为 True 和 False，本流程的表写为 TRUE 和 FALSE；两者都读为逻辑值。
as_logical_columns <- function(x) {
  x[] <- lapply(x, function(v) if (is.character(v) && length(v) && all(v %in% c("True","False","TRUE","FALSE",NA))) v %in% c("True","TRUE") else v)
  x
}

# Compare a new table with the archived one. Rows are sorted by the key columns when keys are given,
# numbers match when all.equal() accepts them (relative difference below 1.5e-8, which absorbs the rounding of
# numbers to 15 significant digits when R writes a table), and text must be identical. Columns in ignore are
# left out.
# 比较新表与存档表。给出键列时先按键列排序；数值以 all.equal() 判断（相对差小于 1.5e-8，可吸收 R 写表时把数字
# 保留 15 位有效数字造成的差别），文字必须完全相同。ignore 中的列不参与比较。
compare_table <- function(name, old_path=file.path(archive, name), new_path=file.path(out, name),
                          keys=NULL, ignore=character(), new=NULL, old=NULL) {
  if (is.null(new)) new <- read.csv(new_path, check.names=FALSE)
  if (is.null(old)) old <- read.csv(old_path, check.names=FALSE)
  new <- as_logical_columns(new); old <- as_logical_columns(old)
  columns <- setdiff(intersect(names(old), names(new)), ignore)
  missing <- setdiff(setdiff(names(old), names(new)), ignore)
  result <- if (nrow(new) != nrow(old)) {
    sprintf("DIFFERS: %d rows, archived %d rows / 不一致：%d 行，存档 %d 行", nrow(new), nrow(old), nrow(new), nrow(old))
  } else {
    if (length(keys)) {
      new <- new[do.call(order, unname(as.list(new[keys]))), , drop=FALSE]
      old <- old[do.call(order, unname(as.list(old[keys]))), , drop=FALSE]
    }
    gap <- 0; fails <- character()
    for (column in columns) {
      a <- new[[column]]; b <- old[[column]]
      if (is.numeric(a) && is.numeric(b)) {
        if (!identical(is.na(a), is.na(b))) { fails <- c(fails, paste(column, "(missing values)")); next }
        if (any(!is.na(a))) gap <- max(gap, abs(a - b), na.rm=TRUE)
        if (!isTRUE(all.equal(a, b))) fails <- c(fails, column)
      } else if (!identical(as.character(a), as.character(b))) fails <- c(fails, paste(column, "(text)"))
    }
    if (!length(fails)) sprintf("matches the archived table (largest difference %.2g) / 与存档表一致（最大差值 %.2g）", gap, gap)
    else sprintf("DIFFERS from the archived table in %s (largest difference %.3g) / 与存档表不一致（%s，最大差值 %.3g）",
                 paste(fails, collapse=", "), gap, paste(fails, collapse=", "), gap)
  }
  if (length(missing)) result <- paste0(result, "; columns not written / 未写出的列: ", paste(missing, collapse=", "))
  cat(sprintf("  %s: %s\n", name, result))
  invisible(result)
}

# Compare a new PNG figure with the archived one: identical files, the same pixels, or the share of pixels that differ.
# 比较新 PNG 图与存档图：文件完全相同、像素相同，或不同像素所占比例。
compare_figure <- function(name, old_path=file.path(archive, name), new_path=file.path(out, name)) {
  same_file <- identical(readBin(new_path, "raw", file.size(new_path)), readBin(old_path, "raw", file.size(old_path)))
  result <- if (same_file) "identical to the archived file / 与存档文件完全相同" else {
    a <- png::readPNG(new_path, native=TRUE); b <- png::readPNG(old_path, native=TRUE)
    if (!identical(dim(a), dim(b))) {
      sprintf("DIFFERS: %s pixels, archived %s pixels / 不一致：%s 像素，存档 %s 像素",
              paste(dim(a), collapse="x"), paste(dim(b), collapse="x"), paste(dim(a), collapse="x"), paste(dim(b), collapse="x"))
    } else {
      share <- mean(unclass(a) != unclass(b))
      if (share == 0) "same pixels as the archived figure / 与存档图像素相同"
      else sprintf("DIFFERS in %.4f%% of pixels / 与存档图有 %.4f%% 的像素不同", 100*share, 100*share)
    }
  }
  cat(sprintf("  %s: %s\n", name, result))
  invisible(result)
}

# ---- Sorting and checks / 排序与核对 ----

# Sort the rows of a table by the key columns in character code order (as the production scripts sort them),
# which does not depend on the language settings of the computer.
# 按键列以字符编码顺序对表的各行排序（与生产脚本的排序相同），不受计算机语言设置的影响。
sort_rows <- function(x, keys) {
  x <- x[do.call(order, c(unname(as.list(x[keys])), method="radix")), , drop=FALSE]
  rownames(x) <- NULL
  x
}

# The character code order of a vector of names. / 一组名称的字符编码顺序。
sorted <- function(x) sort(unique(x), method="radix")

# One text key per row, used to match rows of two tables. / 每行一个文字键，用于匹配两张表的行。
row_key <- function(x, keys) do.call(paste, c(unname(as.list(x[keys])), sep="\r"))

# TRUE when every |a - b| <= atol + rtol |b|; missing values must be missing in both when equal_nan is TRUE.
# 当所有 |a - b| <= atol + rtol |b| 时为 TRUE；equal_nan 为 TRUE 时，缺失值必须在两者中同时缺失。
all_close <- function(a, b, rtol, atol=0, equal_nan=FALSE) {
  both <- is.na(a) & is.na(b)
  near <- !is.na(a) & !is.na(b) & abs(a - b) <= atol + rtol*abs(b)
  length(a) == length(b) && all(near | (equal_nan & both))
}

# ---- Bootstrap samples / 自助抽样 ----
# The production scripts drew the bootstrap samples with the default random number generator of NumPy
# (PCG64 seeded through SeedSequence) and its integers() function. numpy_rng() and bootstrap_draws() repeat
# that generator in base R, so the same seed gives the same samples and the same intervals. Numbers of 32 bits
# are held as doubles and numbers of 128 bits as eight 16-bit parts, lowest first, so that every product stays
# exact.
# 生产脚本以 NumPy 默认随机数生成器（经 SeedSequence 设定种子的 PCG64）及其 integers() 函数进行自助抽样。
# numpy_rng() 与 bootstrap_draws() 在基础 R 中重复该生成器，因此相同种子给出相同的抽样和相同的区间。32 位数以
# 双精度数保存，128 位数以八个 16 位部分（低位在前）保存，使每次乘法都保持精确。

# Multiplication modulo 2^32 and bitwise exclusive or of 32-bit numbers. / 32 位数模 2^32 的乘法与按位异或。
u32_mul <- function(a, b) {
  al <- a %% 65536; ah <- a %/% 65536; bl <- b %% 65536; bh <- b %/% 65536
  (al*bl + ((ah*bl + al*bh) %% 65536)*65536) %% 4294967296
}
u32_xor <- function(a, b) {
  bitwXor(as.integer(a %/% 65536), as.integer(b %/% 65536))*65536 + bitwXor(as.integer(a %% 65536), as.integer(b %% 65536))
}

# Carry, multiplication and addition modulo 2^128; each row of a is one number, b is one number.
# 模 2^128 的进位、乘法与加法；a 的每一行是一个数，b 是一个数。
limb_carry <- function(x) {
  for (k in 1:7) { x[, k+1] <- x[, k+1] + x[, k] %/% 65536; x[, k] <- x[, k] %% 65536 }
  x[, 8] <- x[, 8] %% 65536
  x
}
limb_mul <- function(a, b) {
  r <- matrix(0, nrow(a), 8)
  for (k in 1:8) for (i in 1:k) if (b[k-i+1] != 0) r[, k] <- r[, k] + a[, i]*b[k-i+1]
  limb_carry(r)
}
limb_add <- function(a, b) limb_carry(sweep(a, 2, b, "+"))
limbs_from_u32 <- function(w) as.vector(rbind(w %% 65536, w %/% 65536))

# The generator of numpy.random.default_rng(seed): SeedSequence mixes the seed into four 32-bit words and
# expands them to the initial state and increment of PCG64.
# numpy.random.default_rng(seed) 的生成器：SeedSequence 把种子混合为四个 32 位字，再展开为 PCG64 的初始状态和增量。
numpy_rng <- function(seed) {
  init_a <- 0x43b0d7e5; mult_a <- 0x931e8875; init_b <- 0x8b51f9dd; mult_b <- 0x58f38ded
  entropy <- numeric(); n <- seed
  repeat { entropy <- c(entropy, n %% 4294967296); n <- n %/% 4294967296; if (n == 0) break }
  hash <- init_a
  hashmix <- function(v) {
    v <- u32_xor(v, hash); hash <<- u32_mul(hash, mult_a); v <- u32_mul(v, hash); u32_xor(v, v %/% 65536)
  }
  mix <- function(x, y) { r <- (u32_mul(0xca01f9dd, x) - u32_mul(0x4973f715, y)) %% 4294967296; u32_xor(r, r %/% 65536) }
  pool <- numeric(4)
  for (i in 1:4) pool[i] <- hashmix(if (i <= length(entropy)) entropy[i] else 0)
  for (s in 1:4) for (d in 1:4) if (s != d) pool[d] <- mix(pool[d], hashmix(pool[s]))
  for (s in seq_along(entropy)[-(1:4)]) for (d in 1:4) pool[d] <- mix(pool[d], hashmix(entropy[s]))
  words <- numeric(8)
  for (i in 1:8) {
    v <- u32_xor(pool[(i-1) %% 4 + 1], init_b); init_b <- u32_mul(init_b, mult_b)
    v <- u32_mul(v, init_b); words[i] <- u32_xor(v, v %/% 65536)
  }
  g <- new.env()
  g$mult <- c(0xf645, 0x9fcc, 0xdf64, 0x4385, 0x5da4, 0x1fc6, 0xed05, 0x2360)
  g$inc <- as.vector(limb_carry(matrix(2*limbs_from_u32(words[c(7, 8, 5, 6)]), 1))) + c(1, rep(0, 7))
  state <- limb_add(matrix(g$inc, 1), limbs_from_u32(words[c(3, 4, 1, 2)]))
  g$state <- limb_add(limb_mul(state, g$mult), g$inc)
  g$has32 <- FALSE; g$buf <- 0
  g
}

# The next n states, computed in blocks that double in length (state times multiplier plus increment).
# 接下来的 n 个状态，以长度逐次加倍的块计算（状态乘以乘数再加增量）。
pcg_states <- function(g, n) {
  a <- g$mult; c <- g$inc
  block <- limb_add(limb_mul(g$state, a), c)
  while (nrow(block) < n) {
    block <- rbind(block, limb_add(limb_mul(block, a), c))
    c <- as.vector(limb_add(limb_mul(matrix(c, 1), a), c)); a <- as.vector(limb_mul(matrix(a, 1), a))
  }
  block <- block[seq_len(n), , drop=FALSE]
  g$state <- block[n, , drop=FALSE]
  block
}

# The 32-bit numbers that the generator gives next: each 64-bit output is the exclusive or of the two halves of
# the state, rotated right by its top six bits, and is used as its lower and then its upper 32 bits.
# 生成器接下来给出的 32 位数：每个 64 位输出是状态两半的异或，再按状态最高六位循环右移，先用其低 32 位，再用高 32 位。
next_u32 <- function(g, m) {
  out <- numeric()
  if (m > 0 && g$has32) { out <- g$buf; g$has32 <- FALSE; m <- m - 1 }
  if (m > 0) {
    k <- ceiling(m/2); s <- pcg_states(g, k)
    x <- matrix(bitwXor(as.integer(s[, 1:4]), as.integer(s[, 5:8])), ncol=4)
    turn <- s[, 8] %/% 1024; q <- turn %/% 16; b <- turn %% 16; rows <- seq_len(k)
    y <- matrix(sapply(0:3, function(i) x[cbind(rows, (i + q) %% 4 + 1)]), ncol=4)
    z <- matrix(sapply(1:4, function(i) y[, i] %/% 2^b + (y[, i %% 4 + 1] %% 2^b)*2^(16 - b)), ncol=4)
    lo <- z[, 1] + z[, 2]*65536; hi <- z[, 3] + z[, 4]*65536
    if (m %% 2 == 1) { g$has32 <- TRUE; g$buf <- hi[k] }
    out <- c(out, as.vector(rbind(lo, hi))[seq_len(m)])
  }
  out
}

# Bootstrap samples of n replicates, as rng.integers(n, size=(times, n)) plus one: each row is one sample of
# row positions. Lemire's method maps each 32-bit number to 0, ..., n - 1 and rejects the few that would bias it.
# n 次重复的自助抽样，等同于 rng.integers(n, size=(times, n)) 加一：每行是一次抽样的行位置。Lemire 方法把每个
# 32 位数映射到 0, ..., n - 1，并舍弃少数会造成偏差的数。
bootstrap_draws <- function(g, n, times=10000) {
  stopifnot(n >= 1, n < 2^21)
  if (n == 1) return(matrix(1, times, 1))
  threshold <- (4294967296 - n) %% n
  draws <- numeric(); need <- times*n
  while (need > 0) {
    m <- next_u32(g, need)*n
    keep <- m %% 4294967296 >= threshold
    draws <- c(draws, (m %/% 4294967296)[keep]); need <- need - sum(keep)
  }
  matrix(draws + 1, times, n, byrow=TRUE)
}

# The 2.5% and 97.5% quantiles of the medians of the bootstrap samples (type 7, as numpy.quantile).
# 各次自助抽样中位数的 2.5% 与 97.5% 分位数（第 7 类，与 numpy.quantile 相同）。
bootstrap_interval <- function(values, draws) {
  boot <- apply(matrix(values[draws], nrow(draws)), 1, median)
  unname(quantile(boot, c(.025, .975), type=7))
}
