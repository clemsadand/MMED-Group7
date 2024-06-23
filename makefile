# Define variables
REPO_URL = https://github.com/clemsadand/MMED-Group7.git
FILENAME = .

# Define targets
.PHONY: push

# Push the file to GitHub
push:
	git add $(FILENAME)
	git commit -m "This version is based on the last ODE equations(see Overleaf)"
	git push origin main
