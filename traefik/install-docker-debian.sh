for pkg in docker.io docker-doc docker-compose podman-docker containerd runc; do sudo apt-get remove $pkg; done

# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/debian/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/debian \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update

sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

sudo echo "GNNEPCSAFT_SECRET_KEY=x5F1GrWg1gEc1foE3CpPTi7b6h3CfcFbXmdJ00Us5kNYmo7TcG
GNNEPCSAFT_DOMAIN_NAME=https://*.wildsonbbl.com
GNNEPCSAFT_ALLOWED_HOSTS=.wildsonbbl.com, .localhost, .0.0.0.0
GOOGLE_API_KEY=notvalidapikey
GNNEPCSAFT_DEBUG=false
DJANGO_SETTINGS_MODULE=webapp.settings
GNNEPCSAFT_LOG_LEVEL=WARNING
GNNEPCSAFT_PLATFORM=webapp" >> .traefik-env

sudo docker run -d -p 80:8000 --env-file .traefik-env --name webapp  wildsonbbl/gnnepcsaft