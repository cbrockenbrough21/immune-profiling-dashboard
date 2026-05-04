.PHONY: setup generate pipeline dashboard

setup:
	pip install -r requirements.txt

generate:
	python generate_data.py

pipeline:
	python load_data.py
	python run_analysis.py

dashboard:
	streamlit run dashboard.py