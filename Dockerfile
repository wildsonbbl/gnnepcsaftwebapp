# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3.9.18-slim-bookworm
EXPOSE 8000

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install curl unzip gcc g++ -y

# Install pip requirements
RUN python -m pip install -U pip wheel
COPY requirements.txt .
RUN python -m pip install -r requirements.txt
RUN curl -O -L https://gitlab.com/libeigen/eigen/-/archive/master/eigen-master.zip
RUN curl -O -L https://github.com/zmeri/PC-SAFT/archive/refs/tags/v1.4.1.zip
RUN unzip -q eigen-master.zip
RUN unzip -q v1.4.1.zip
RUN cp -rf eigen-master/. PC-SAFT-1.4.1/externals/eigen
RUN python -m pip install ./PC-SAFT-1.4.1 --no-build-isolation
RUN rm -rf PC-SAFT* eigen-master* v1.4.1.zip

RUN python manage.py collectstatic --no-input
WORKDIR /app
COPY . /app

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
CMD ["gunicorn", "--bind", "0.0.0.0:8000", "webapp.wsgi"]
