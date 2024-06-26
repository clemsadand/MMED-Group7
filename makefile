# Define variables
REPO_URL = https://github.com/clemsadand/MMED-Group7.git
FILENAME = code20240626.R

# Define targets
.PHONY: push

# Push the file to GitHub
push:
	git add .
	git commit -m "With intervention"
	git push origin main
