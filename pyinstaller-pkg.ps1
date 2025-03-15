pyinstaller --distpath ./apppkg/dist --workpath ./apppkg/build `
  -D -n gnnepcsaftwebapp --additional-hooks-dir="./hooks" --noconfirm --clean `
  --add-data="./icons.json:." `
  --add-data="./mydatabase:." `
  --add-data="./gnnmodel/templates:./gnnmodel/templates" `
  --add-data="./productionfiles:./productionfiles" `
  --exclude-module="polars" `
  --copy-metadata django-bootstrap-v5 `
  manage.py