import re

with open("api.md") as f:
    text = f.read()

# Pattern: match anchor lines like
# <a id="circuits.Circuits.__init__"></a>
pattern = r"(<a id=\".+?\"></a>)"

# Replace with itself + horizontal line
text = re.sub(pattern, r"\1\n\n---\n", text)

with open("api.md", "w") as f:
    f.write(text)
