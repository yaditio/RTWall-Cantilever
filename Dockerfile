# Use official Python 3.10 slim base image
FROM python:3.10-slim

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1

# Set working directory
WORKDIR /app

# Install system dependencies for OpenSees, sectionproperties, and CairoSVG
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3-dev \
    libffi-dev \
    gfortran \
    libgl1 \
    libglu1-mesa \
    libgomp1 \
    libcairo2 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf2.0-0 \
    shared-mime-info \
    && rm -rf /var/lib/apt/lists/*

# Copy dependencies list
COPY requirements.txt .

# Install Python packages
RUN pip install --no-cache-dir -U pip && \
    pip install --no-cache-dir -r requirements.txt

# Copy application source code
COPY . .

# Expose default Streamlit port
EXPOSE 8501

# Command to run Streamlit on container startup
CMD ["streamlit", "run", "App.py", "--server.port=8501", "--server.address=0.0.0.0"]
