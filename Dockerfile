# Use R base image
FROM rocker/r-ver:4.5.2
# Install system dependencies including JAGS
RUN apt-get update && apt-get install -y --no-install-recommends jags && rm -rf /var/lib/apt/lists/*
# Set working directory
WORKDIR /app
# Copy all project files
COPY . .

# Set environment variables so R can find project library from any working directory
ENV R_PROFILE_USER=/app/.Rprofile

# Install renv and restore dependencies
RUN R -e "install.packages('renv'); renv::consent(provided = TRUE); renv::restore(prompt = FALSE)"