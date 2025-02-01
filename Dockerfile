# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3.10-bookworm
EXPOSE 8000

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install curl unzip gcc g++ libxrender1 libxext6 -y

# Install pip requirements
RUN python -m pip install -U pip wheel
COPY requirements.txt .
RUN python -m pip install -r requirements.txt

WORKDIR /app
COPY . /app
COPY daemon.json /etc/docker/daemon.json
RUN python manage.py collectstatic --no-input
RUN python manage.py migrate --no-input

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
CMD ["gunicorn", "--bind", "0.0.0.0:8000", "webapp.wsgi"]
