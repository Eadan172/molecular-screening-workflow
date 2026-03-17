import pickle

# 加载模型
model_path = './models/ultimate_ensemble_qsar_model.pkl'
with open(model_path, 'rb') as f:
    model_data = pickle.load(f)

print(f"模型数据类型: {type(model_data)}")
print(f"模型数据长度: {len(model_data)}")
for i, item in enumerate(model_data):
    print(f"第{i}项类型: {type(item)}")
