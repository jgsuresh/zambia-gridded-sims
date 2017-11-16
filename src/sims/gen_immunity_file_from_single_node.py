# Inputs:
#   --demographics file (e.g. generated in example_josh)
#   --single-node immunity file (e.g. "Sinamalima_1_node_immune_init_p1_33_p2_117.json")

# Outputs:
#   --multi-node immunity file where each node has identical immunity values.

import json

demo_fp = "./inputs/Demographics/demo_test.json"
imm_1node_fp = "./inputs/Immunity/Sinamalima_1_node_immune_init_p1_33_p2_117.json"
imm_multinode_fp = "./inputs/Immunity/immune_test_p1_33_p2_117.json"

# Open demographics file
with open(demo_fp,'r') as f:
    demo_dict = json.load(f)

# Get NodeID list
node_ids = []
for node in demo_dict['Nodes']:
    node_ids.append(node['NodeID'])

# Open single-node immunity file
with open(imm_1node_fp,'r') as f:
    imm_dict = json.load(f)

# Edit node list in this dictionary
imm_node_list = imm_dict['Nodes']
for node_id in node_ids:
    imm_node_list.append({u'NodeID':node_id})

del imm_node_list[0] # Remove the original node that was in the single-node immunity file

# Edit node metadata to reflect new number of nodes:
imm_dict['Metadata']['NodeCount'] = len(imm_node_list)

# Dump as new json file
with open(imm_multinode_fp,'w') as f:
    json.dump(imm_dict,f,indent=4)