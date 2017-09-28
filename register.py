try:
    long_description = pypandoc.convert('README.md', 'rst')
    long_description = long_description.replace("\r","") # YOU  NEED THIS LINE
except OSError:
    print("Pandoc not found. Long_description conversion failure.")
    import io
    # pandoc is not installed, fallback to using raw contents
    with io.open('README.md', encoding="utf-8") as f:
        long_description = f.read()