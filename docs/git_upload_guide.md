# Git上传指引

本文档提供详细的Git仓库初始化和代码上传步骤。

## 前置准备

### 1. 安装Git

**Windows**:
- 下载并安装: https://git-scm.com/download/win
- 或使用 `winget install Git.Git`

**Mac**:
```bash
brew install git
```

**Linux**:
```bash
sudo apt-get install git  # Ubuntu/Debian
sudo yum install git      # CentOS/RHEL
```

### 2. 配置Git

```bash
# 设置用户名和邮箱
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"

# 验证配置
git config --list
```

### 3. 创建GitHub仓库

1. 登录 GitHub (https://github.com)
2. 点击右上角 "+" → "New repository"
3. 填写仓库信息:
   - Repository name: `molecular-screening-workflow`
   - Description: `基于RNN的分子生成与多阶段虚拟筛选工作流`
   - 选择 Public 或 Private
   - **不要**勾选 "Initialize this repository with a README"（我们已有README）
4. 点击 "Create repository"

## 上传步骤

### 方法1: 从现有目录初始化（推荐）

```bash
# 1. 进入项目目录
cd molecular-screening-workflow

# 2. 初始化Git仓库
git init

# 3. 添加所有文件到暂存区
git add .

# 4. 创建首次提交
git commit -m "Initial commit: 分子筛选工作流完整代码

- 添加统一入口程序 run_workflow.py
- 添加模块化源代码 (src/)
- 添加配置文件 (config/)
- 添加Demo演示数据 (demo/)
- 添加详细文档 (docs/)
- 添加README和依赖列表"

# 5. 添加远程仓库
git remote add origin https://github.com/yourusername/molecular-screening-workflow.git

# 6. 推送到GitHub
git branch -M main
git push -u origin main
```

### 方法2: 使用GitHub CLI（更简单）

```bash
# 1. 安装GitHub CLI
# Windows: winget install GitHub.cli
# Mac: brew install gh
# Linux: sudo apt install gh

# 2. 登录GitHub
gh auth login

# 3. 进入项目目录
cd molecular-screening-workflow

# 4. 初始化并创建仓库
git init
git add .
git commit -m "Initial commit: 分子筛选工作流完整代码"
gh repo create molecular-screening-workflow --public --source=. --push
```

## 完整命令示例

```bash
#!/bin/bash
# Git上传完整脚本

# 进入项目目录
cd /path/to/molecular-screening-workflow

# 初始化
git init

# 添加文件
git add .

# 查看状态
git status

# 提交
git commit -m "Initial commit: 分子筛选工作流完整代码

功能模块:
- RNN分子生成
- QSAR活性预测
- ADMET性质过滤
- 分子对接
- SA Score评估
- 结果可视化

文档:
- README.md (项目说明)
- docs/ (详细文档)
- demo/ (演示数据)"

# 添加远程仓库（替换为你的用户名）
git remote add origin https://github.com/YOUR_USERNAME/molecular-screening-workflow.git

# 设置主分支
git branch -M main

# 推送
git push -u origin main

echo "✅ 代码已成功上传到GitHub！"
```

## 后续更新

当需要更新代码时：

```bash
# 1. 查看修改
git status

# 2. 添加修改的文件
git add .

# 3. 提交修改
git commit -m "描述你的修改"

# 4. 推送到GitHub
git push
```

## 分支管理

### 创建新分支

```bash
# 创建并切换到新分支
git checkout -b feature/new-feature

# 在新分支上工作
git add .
git commit -m "Add new feature"

# 推送新分支
git push origin feature/new-feature
```

### 合并分支

```bash
# 切换回主分支
git checkout main

# 合并新分支
git merge feature/new-feature

# 推送合并结果
git push
```

## .gitignore说明

项目已配置 `.gitignore` 文件，以下内容不会被上传：

- `__pycache__/` - Python缓存
- `*.pyc` - 编译文件
- `venv/` - 虚拟环境
- `logs/` - 日志文件
- `results/**/*.csv` - 结果文件
- `models/*.pkl` - 大型模型文件
- `.DS_Store` - Mac系统文件

## 大文件处理

如果模型文件过大（>100MB），建议使用Git LFS：

```bash
# 1. 安装Git LFS
git lfs install

# 2. 跟踪大文件
git lfs track "*.pkl"
git lfs track "*.keras"
git lfs track "*.h5"

# 3. 提交.gitattributes
git add .gitattributes
git commit -m "Configure Git LFS"

# 4. 正常提交大文件
git add models/
git commit -m "Add model files"
git push
```

## 常见问题

### Q1: 推送失败，提示"fatal: 'origin' already exists"

**解决方案**:
```bash
git remote remove origin
git remote add origin https://github.com/yourusername/molecular-screening-workflow.git
```

### Q2: 推送被拒绝，提示"Updates were rejected"

**解决方案**:
```bash
# 强制推送（谨慎使用，会覆盖远程内容）
git push -f origin main

# 或先拉取再推送
git pull origin main --allow-unrelated-histories
git push origin main
```

### Q3: 如何删除已上传的敏感文件？

**解决方案**:
```bash
# 从Git历史中删除文件
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch path/to/sensitive-file" \
  --prune-empty --tag-name-filter cat -- --all

# 强制推送
git push origin --force --all
```

### Q4: 如何撤销最近的提交？

**解决方案**:
```bash
# 撤销最近一次提交，保留修改
git reset --soft HEAD~1

# 撤销最近一次提交，丢弃修改
git reset --hard HEAD~1
```

## GitHub仓库设置建议

上传完成后，建议进行以下设置：

1. **添加Topics标签**:
   - `drug-discovery`
   - `molecular-docking`
   - `machine-learning`
   - `cheminformatics`
   - `qsar`

2. **创建Release**:
   - 点击 "Releases" → "Draft a new release"
   - 标签: `v1.0.0`
   - 标题: `Initial Release`
   - 描述: 首次发布，包含完整工作流

3. **启用GitHub Pages**（可选）:
   - Settings → Pages
   - 选择分支和目录
   - 用于托管文档网站

4. **添加徽章**:
   在README.md顶部添加：
   ```markdown
   [![GitHub stars](https://img.shields.io/github/stars/yourusername/molecular-screening-workflow.svg)](https://github.com/yourusername/molecular-screening-workflow/stargazers)
   [![GitHub forks](https://img.shields.io/github/forks/yourusername/molecular-screening-workflow.svg)](https://github.com/yourusername/molecular-screening-workflow/network)
   [![GitHub issues](https://img.shields.io/github/issues/yourusername/molecular-screening-workflow.svg)](https://github.com/yourusername/molecular-screening-workflow/issues)
   ```

## 验证上传

上传完成后，访问你的仓库页面验证：

1. ✅ 所有文件都已上传
2. ✅ README.md正确显示
3. ✅ 目录结构清晰
4. ✅ .gitignore生效（临时文件未上传）

## 获取仓库链接

上传成功后，你的仓库地址为：
```
https://github.com/yourusername/molecular-screening-workflow
```

其他人可以通过以下命令克隆仓库：
```bash
git clone https://github.com/yourusername/molecular-screening-workflow.git
```

---

**恭喜！你的分子筛选工作流已成功上传到GitHub！** 🎉
