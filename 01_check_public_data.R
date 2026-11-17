# =====================================================================================================
# Step 01. Check the stored copy of the public IATTC longline file (Section 2.1)
# 第 01 步  核对已保存的 IATTC 公开延绳钓文件（第 2.1 节）
#
# The analysis uses one public file, PublicLLTunaBillfishNum.csv: the monthly public domain records of
# retained catch in numbers and numbers of hooks by flag and 5-degree cell, from the IATTC archive
# PublicLLTunaBillfish.zip (release "March 2026", downloaded on 2026-11-06 into data/raw/). This step checks
# that the stored CSV is the same as the copy inside the stored zip archive. IATTC replaces the file at each
# release, so the walkthrough does not download it again. The raw data are only read, never changed.
# 本分析只使用一个公开文件 PublicLLTunaBillfishNum.csv：按船旗和 5° 网格统计的月度延绳钓公开记录，包括留存
# 尾数和钩数，来自 IATTC 压缩包 PublicLLTunaBillfish.zip（"March 2026" 版本，2026-11-06 下载到 data/raw/）。
# 本步骤检查保存的 CSV 与保存的压缩包内的文件相同。IATTC 每次发布都会替换该文件，因此本流程不重新下载。原始数据
# 只读取、不修改。
#
# Reads / 读取:  config/public_sources.json, config/directed_fishery.json, data/raw/manifest.json,
#                data/raw/PublicLLTunaBillfish.zip, data/raw/longline/PublicLLTunaBillfishNum.csv
# Writes / 写出: nothing / 无
#
# Run from the project folder; the first line makes R use one core and read the Chinese folder names.
# 在项目文件夹中运行；第一行让 R 只用一个核心，并能读取中文文件夹名。
#   export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
#   Rscript analysis_walkthrough/01_check_public_data.R
# =====================================================================================================

cat("Step 01: checking the stored public file / 第 01 步：核对保存的公开文件\n")

# Read the list of sources, the analysis settings and the record written when the data were downloaded.
# 读取数据来源清单、分析设置和下载数据时写下的记录。
sources <- jsonlite::fromJSON("config/public_sources.json")$sources
cfg <- jsonlite::fromJSON("config/directed_fishery.json")
manifest <- jsonlite::fromJSON("data/raw/manifest.json")
source_row <- sources[sources$filename == "PublicLLTunaBillfish.zip", ][1, ]
cat("Source / 数据来源:", source_row$url, "\n")
cat("Downloaded / 下载日期:", manifest$access_date_local, "| release / 版本:", source_row$release_on_landing_page, "\n")

# Check that the stored CSV is the same as the copy inside the stored zip archive, byte for byte.
# 逐字节检查保存的 CSV 与保存的压缩包内的文件相同。
zip_path <- file.path("data/raw", source_row$filename)
csv_path <- cfg$num
size <- file.size(csv_path)
connection <- unz(zip_path, basename(csv_path), open="rb")
inside <- readBin(connection, "raw", size + 1)
close(connection)
same <- identical(inside, readBin(csv_path, "raw", size + 1))
cat(sprintf("  %s: %s\n", basename(csv_path), if (same) "same as the copy in the zip archive / 与压缩包内的文件相同"
            else "DIFFERS from the copy in the zip archive / 与压缩包内的文件不同"))

# Count the lines: a line ends with a line feed, a carriage return, or both; the first line is the header.
# 统计行数：每行以换行符、回车符或两者结尾；第一行为表头。
lf <- inside == as.raw(10); cr <- inside == as.raw(13)
breaks <- sum(lf) + sum(cr & !c(lf[-1], FALSE))
lines <- breaks + (length(inside) > 0 && !(lf[length(inside)] || cr[length(inside)]))
cat("Rows in the file / 文件行数:", lines - 1, "\n")
