# 公开发布前执行清单

## 必须完成（未完成不要发布v1.0.0）

- [ ] 填写发现参考基因组和GFF3的公共accession、版本、URL与checksum。
- [ ] 确认发现样本原始reads是否为原始下机数据或公司已清洗数据，并补齐清洗参数。
- [ ] 填写60个发现样本和60个外部样本的BioSample/Run accession。
- [ ] 补齐采样时间、地点、文库构建、平台、读长、伦理批准号。
- [ ] 说明LG11与CM056831.1之间的坐标/等位基因匹配方法并发布crosswalk。
- [ ] 恢复PLINK、ADMIXTURE完整版本和命令。
- [ ] 用VEP或SnpEff复核编码后果。
- [ ] 把最终绘图脚本改为只读取仓库相对路径；运行全部图表的clean-environment测试。
- [ ] 确认主文和补充材料采用60尾性腺解剖真值：SNP为54/60（90.0%）；严格k-mer为56/58（96.6%），覆盖58/60；BY1/BY5保留为不确定。
- [ ] 选择并写入代码和数据许可证。
- [ ] 删除所有`AUTHOR_INPUT_REQUIRED`。

## 数据库提交顺序

1. 建立BioProject。
2. 为每尾鱼建立BioSample，不把文库重复误作生物学个体。
3. 准备FASTQ md5和SRA/GSA样本表。
4. 上传发现与外部验证reads，取得Run accession。
5. 将accession回填到仓库样本表和论文Data Availability。
6. 将代码推送到私有GitHub，运行本仓库审计。
7. 创建`v1.0.0` release，并由Zenodo归档获得DOI。
8. 将DOI写入论文、README和CITATION.cff。

## 一键审计

```bash
bash workflow/run_source_data_checks.sh
```

当前release candidate预期返回退出码2，因为仍存在作者必须补充的字段；只有退出码0才可视为正式公开版本。
