# Define variables
REPO_URL = https://github.com/clemsadand/MMED-Group7.git
FILENAME = TheCode.R

# Define targets
.PHONY: push

# Push the file to GitHub
push:
	git add .
	git commit -m "Add for data simulation"
	git push origin main
