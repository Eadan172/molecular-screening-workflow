# 分子设计流程

下载代码后，在项目目录执行一条命令创建 `.venv`。把分子要求写进 `input/request.txt` 即可运行。

填写 `.env` 里的 `LLM_API_KEY` 后，程序会调用兼容 OpenAI 接口的大模型，完成要求分析、任务拆解补充、调研说明、计算数据解释和报告整理。不填密钥时，同样的命令会跳过大模型，直接做规则解析、公开数据库检索和本地计算，并写出规则报告。

## 一条命令准备环境

Linux 或 macOS：

```bash
bash setup.sh
```

Windows：

```bat
setup.bat
```

脚本会在当前项目目录创建 `.venv`，安装 RDKit 和 requests，并在没有 `.env` 时从 `.env.example` 复制一份。终端里直接回车可以跳过密钥。

然后编辑两个文件：

- `.env`：只在需要大模型时填写 `LLM_API_KEY`、`LLM_BASE_URL`、`LLM_MODEL`
- `input/request.txt`：靶点、种属、理化性质、活性要求和可选的计算软件

运行：

```bash
bash run.sh
```

Windows 使用 `run.bat`。`run.sh` / `run.bat` 发现还没有 `.venv` 时会先安装环境。

只跑计算、即使已经写了密钥也不调用模型：

```bash
bash run.sh --calc-only
```

指定别的需求文件：

```bash
bash run.sh input/my_request.txt
```

结果写到 `results/runs/<时间戳>/`，其中 `report.md` 是报告，`properties.csv` 是全部分子性质。

## 输入文件写什么

用方括号分段，用「键: 值」写要求。`#` 开头是注释。示例见 `input/request.txt`。

```text
[基本信息]
靶点: AKT1
种属: Homo sapiens
适应症: 肿瘤
目标: 寻找可口服的小分子抑制剂

[理化性质]
分子量: <= 500
LogP: 1-3
TPSA: <= 140
氢键供体: <= 5
氢键受体: <= 10
可旋转键: <= 10
QED: >= 0.3

[活性与成药性]
IC50: < 100 nM
hERG: 避免
AMES: 阴性

[分子]
CC(=O)Oc1ccccc1C(=O)O
分子文件: data/candidates.smi

[计算软件]
name: vina
type: local
path: /usr/local/bin/vina
receptor: pdb_files_AKT1/4gv1.pdb
args: --receptor {receptor} --ligand {ligand} --out {output}

name: predict_api
type: api
url: http://127.0.0.1:8000/predict
method: POST
header: Authorization: Bearer YOUR_TOOL_TOKEN
```

没有 `[分子]` 时，程序仍然分析要求和检索靶点，但不会做性质计算。

数值可以写成 `<= 500`、`< 100 nM` 或 `1-3`。只写一个数字时，QED 和 Fsp3 按大于等于处理，其余按小于等于处理。

内置计算使用 RDKit：分子量、LogP、TPSA、氢键供体/受体、可旋转键、Fsp3、QED、Lipinski、Veber，以及能用时的 SA Score。Lipinski 和 Veber 只作为参考列。结构警示只是子结构提示，不能代替 hERG、AMES 或实验。

输入里写了 IC50、Ki 或结合能，但没有对应数据时，分子会标成「数据不足」，不会被假装成已经通过。外部程序或接口如果返回带 `smiles` 列，以及 `ic50` 或 `binding_affinity` 的 CSV/JSON，程序会按 SMILES 对齐后再判断。

## 本地软件和 API

`[计算软件]` 里每一段以 `name:` 开始。

本地程序：

- `type: local`
- `path`: 可执行文件的绝对路径、相对项目的路径，或 `PATH` 里的命令名
- `args`: 参数。先按空格拆分，再替换占位符，因此路径里可以有空格
- `receptor`: 受体文件
- `timeout`: 秒，默认 300

占位符：

| 占位符 | 含义 |
| --- | --- |
| `{ligand}` | 由输入 SMILES 生成的 SDF |
| `{smiles}` | 每行一个 SMILES 的文本 |
| `{receptor}` | 受体文件 |
| `{output}` | 建议的结果文件 `result.out` |
| `{outdir}` | 该工具的工作目录 |

接口：

- `type: api`
- `url`: 以 `http://` 或 `https://` 开头
- `method`: `POST` 或 `GET`
- `header`: 一行请求头，例如 `Authorization: Bearer ...`

程序向接口提交靶点、种属、阈值和 SMILES。请求头只用于这次调用，不会写入报告。路径不存在或接口失败时，这一项记为跳过或失败，其他计算继续。

对接软件如果要求 PDBQT，请把转换和对接写进你自己的命令或脚本，再把脚本路径填到 `path`。

## 有密钥和没有密钥

没有 `LLM_API_KEY`，或密钥仍是 `your_api_key` 这类占位符：

1. 用规则读取靶点、种属、阈值、分子和软件
2. 用 ChEMBL 公开接口检索靶点；网络失败只跳过检索
3. 运行 RDKit 和输入文件里声明的本地程序或接口
4. 按阈值统计通过、未通过和数据不足
5. 写出模板报告

填写密钥后，额外调用三次模型：理解需求、根据数据库摘要写调研说明、解释计算结果。模型不能新增要执行的命令。模型调用失败时，该步退回规则结果，计算仍然保留。

`.env` 示例：

```text
LLM_API_KEY=sk-你的密钥
LLM_BASE_URL=https://api.deepseek.com/v1
LLM_MODEL=deepseek-chat
```

`LLM_BASE_URL` 需要包含 `/v1`。也支持 OpenAI、通义千问兼容模式和 Moonshot，具体地址写在 `.env.example`。

密钥放在 `.env`，这个文件不会被 Git 跟踪。不要把密钥写进 `input/request.txt`。工具自己的请求头可以写在计算软件段，它不会进入报告。

## 报告里有什么

`report.md` 按下面的顺序写：

1. 要求分析
2. 任务拆解
3. 调研
4. 计算
5. 计算数据分析
6. 报告整理

同目录还有 `properties.csv`、`spec.json`、`research.json`、`calculations.json` 和 `run.log`。

## 测试

```bash
python -m unittest discover -s tests -v
```

装好 `.venv` 后，用 `.venv/bin/python -m unittest discover -s tests -v` 可以同时跑 RDKit 计算断言。

## 旧的固定筛选脚本

`run_workflow.py` 和 `src/` 里的 RNN、QSAR、ADMET 脚本是原先的固定流程。当前默认入口是 `run.py`。旧脚本中还有未处理的合并冲突，不作为现在的使用方式。若要训练那些模型，另装 `requirements-ml.txt`；其中的 TensorFlow 和 PyTorch 不在一键安装里。

`docs/` 里的架构说明描述的是旧流程。
