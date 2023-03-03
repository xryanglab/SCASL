from scasl import get_config, SCASL
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-y', '--yaml', help='configuration YAML file', type=str, required=True)
cfg_path = parser.parse_args().yaml
cfg = get_config(cfg_path)

scasl = SCASL(cfg)
scasl.fit()