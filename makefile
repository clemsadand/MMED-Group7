# Define variables
REPO_URL = https://github.com/clemsadand/MMED-Group7.git
FILENAME = code220623.R

# Define targets
.PHONY: push

# Push the file to GitHub
push:
	git add .
	git commit -m "Update code220623.R"
	git push origin main
