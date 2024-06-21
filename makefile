# Define variables
REPO_URL = https://github.com/clemsadand/MMED-Group7.git
FILENAME = .

# Define targets
.PHONY: push

# Push the file to GitHub
push:
	git add $(FILENAME)
	git commit -m "Add $(FILENAME)"
	git push origin main
