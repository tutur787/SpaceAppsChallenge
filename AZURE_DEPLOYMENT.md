# Azure Deployment Guide

This guide explains how to deploy the Impactor-2025 app to Azure with both Next.js (background) and Streamlit running together.

## Prerequisites

- Azure account
- Azure CLI installed (`az` command)
- Git repository pushed to GitHub/GitLab/Azure DevOps

## Deployment Steps

### 1. Create Azure Web App

```bash
# Login to Azure
az login

# Create a resource group
az group create --name impactor-rg --location eastus

# Create an App Service plan (Linux, B1 tier minimum for Node.js + Python)
az appservice plan create \
  --name impactor-plan \
  --resource-group impactor-rg \
  --sku B2 \
  --is-linux

# Create the web app with Python 3.9 runtime
az webapp create \
  --resource-group impactor-rg \
  --plan impactor-plan \
  --name impactor-2025 \
  --runtime "PYTHON:3.9"
```

### 2. Configure Web App Settings

```bash
# Set the startup command to run your script
az webapp config set \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --startup-file "start.sh"

# Add Node.js to the Python container
az webapp config appsettings set \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --settings \
    NODE_VERSION="18" \
    NPM_VERSION="9" \
    SCM_DO_BUILD_DURING_DEPLOYMENT="true" \
    ENABLE_ORYX_BUILD="true"

# Configure port (Azure expects 8000 by default)
az webapp config appsettings set \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --settings PORT="8000"
```

### 3. Deploy from Git

**Option A: Deploy from GitHub**

```bash
# Configure GitHub deployment
az webapp deployment source config \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --repo-url https://github.com/YOUR_USERNAME/SpaceAppsChallenge \
  --branch main \
  --manual-integration
```

**Option B: Deploy from local Git**

```bash
# Enable local git deployment
az webapp deployment source config-local-git \
  --resource-group impactor-rg \
  --name impactor-2025

# Get deployment credentials
az webapp deployment list-publishing-credentials \
  --resource-group impactor-rg \
  --name impactor-2025

# Add Azure remote and push
git remote add azure https://<deployment-username>@impactor-2025.scm.azurewebsites.net/impactor-2025.git
git push azure main
```

**Option C: Deploy via ZIP**

```bash
# Create a ZIP of your project
zip -r app.zip . -x "*.git*" "venv/*" "my-app/node_modules/*" "my-app/.next/*"

# Deploy the ZIP
az webapp deployment source config-zip \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --src app.zip
```

### 4. Monitor Deployment

```bash
# Stream deployment logs
az webapp log tail \
  --resource-group impactor-rg \
  --name impactor-2025

# View application logs
az webapp log tail \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --provider application
```

### 5. Access Your App

Your app will be available at:
```
https://impactor-2025.azurewebsites.net
```

## File Structure for Azure

```
SpaceAppsChallenge/
├── start.sh                    # Azure startup script (runs both Next.js & Streamlit)
├── requirements.txt            # Python dependencies
├── main.py                     # Streamlit app entry point
├── src/
│   ├── background.py           # Background component
│   └── ...
├── my-app/                     # Next.js app
│   ├── package.json
│   ├── next.config.ts
│   └── ...
└── .gitignore                  # Excludes node_modules, .next, venv, etc.
```

## Important Azure Configuration

### startup script (`start.sh`)
- Installs Python dependencies from `requirements.txt`
- Installs Node.js dependencies with `npm install --legacy-peer-deps`
- Starts Next.js on port 3000 in background
- Starts Streamlit on Azure's `$PORT` (defaults to 8000)

### Environment Variables (if needed)
Add via Azure Portal or CLI:

```bash
az webapp config appsettings set \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --settings \
    NASA_API_KEY="your_key_here" \
    MAPBOX_API_KEY="your_key_here"
```

## Troubleshooting

### Logs Not Showing
```bash
# Enable application logging
az webapp log config \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --application-logging filesystem \
  --level verbose
```

### Next.js Not Starting
- Check that Node.js 18+ is available: `az webapp config show`
- Verify `npm install` succeeds in logs
- Ensure `start.sh` has execute permissions (should be committed with `chmod +x`)

### Streamlit Not Accessible
- Verify `PORT` environment variable is set to 8000
- Check that Streamlit starts with `--server.address=0.0.0.0`
- Review logs for port binding errors

### Background Not Loading
- Confirm Next.js is running on `http://localhost:3000` (check logs)
- Ensure iframe points to `localhost:3000` (not external URL)
- Both services must run on same Azure instance

## Scaling

For production, upgrade to a higher tier:

```bash
az appservice plan update \
  --name impactor-plan \
  --resource-group impactor-rg \
  --sku S1  # Standard tier
```

## Cost Optimization

- **Development**: B1 tier (~$13/month)
- **Production**: S1 tier (~$70/month) with auto-scaling
- **Stop when not in use**:
  ```bash
  az webapp stop --resource-group impactor-rg --name impactor-2025
  ```

## Security

### Enable HTTPS only
```bash
az webapp update \
  --resource-group impactor-rg \
  --name impactor-2025 \
  --https-only true
```

### Add custom domain (optional)
```bash
az webapp config hostname add \
  --resource-group impactor-rg \
  --webapp-name impactor-2025 \
  --hostname www.yourdomain.com
```

## Additional Resources

- [Azure App Service Docs](https://docs.microsoft.com/azure/app-service/)
- [Deploy Python to Azure](https://docs.microsoft.com/azure/app-service/quickstart-python)
- [Streamlit Deployment Guide](https://docs.streamlit.io/streamlit-community-cloud/get-started/deploy-an-app)
