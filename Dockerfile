FROM python:3.12-slim

WORKDIR /opt/insitugame
COPY . .

RUN python -m pip install --no-cache-dir .

ENTRYPOINT ["insitugame"]
