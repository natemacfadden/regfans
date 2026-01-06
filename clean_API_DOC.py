import re

with open("API_DOC.md") as f:
    text = f.read()

# Pattern: match anchor lines like
# <a id="circuits.Circuits.__init__"></a>
pattern = r"(<a id=\".+?\"></a>)"

# Replace with itself + horizontal line
text = re.sub(pattern, r"\1\n\n---\n", text)

with open("API_DOC.md", "w") as f:
    f.write(text)
