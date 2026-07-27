FROM python:3.11-slim

WORKDIR /app

# No apt layer. Every runtime dependency ships a manylinux wheel, and the two
# fonts the social card needs are committed under static/fonts precisely so
# this image does not have to install any.

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

EXPOSE 8080
ENV PORT=8080 PYTHONUNBUFFERED=1

# gunicorn, not the Flask development server. --preload parses the 7 MB
# proteome once in the parent so both workers share it copy-on-write instead of
# holding a private 11 MB copy each.
CMD ["sh", "-c", "gunicorn --bind 0.0.0.0:$PORT --workers 2 --threads 4 --timeout 60 --preload app:app"]
