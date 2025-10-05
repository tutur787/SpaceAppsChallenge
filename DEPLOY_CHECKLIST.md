# Azure Deployment Quick Checklist ✅

## Before You Deploy

- [ ] Both Next.js and Streamlit work locally
  ```bash
  # Terminal 1
  cd my-app && npm run dev

  # Terminal 2
  streamlit run main.py
  ```

- [ ] Test the combined startup script locally
  ```bash
  ./start.sh
  ```

- [ ] Commit all changes to git
  ```bash
  git add .
  git commit -m "Ready for Azure deployment"
  git push origin main
  ```

## Azure Setup (One-Time)

1. **Create Azure resources:**
   ```bash
   az login
   az group create --name impactor-rg --location eastus
   az appservice plan create --name impactor-plan --resource-group impactor-rg --sku B2 --is-linux
   az webapp create --resource-group impactor-rg --plan impactor-plan --name impactor-2025 --runtime "PYTHON:3.9"
   ```

2. **Configure the app:**
   ```bash
   az webapp config set --resource-group impactor-rg --name impactor-2025 --startup-file "start.sh"

   az webapp config appsettings set --resource-group impactor-rg --name impactor-2025 \
     --settings NODE_VERSION="18" PORT="8000" ENABLE_ORYX_BUILD="true"
   ```

3. **Deploy your code:**

   **Option 1 - From GitHub (recommended):**
   ```bash
   az webapp deployment source config --resource-group impactor-rg --name impactor-2025 \
     --repo-url https://github.com/YOUR_USERNAME/SpaceAppsChallenge --branch main --manual-integration
   ```

   **Option 2 - From local:**
   ```bash
   az webapp deployment source config-local-git --resource-group impactor-rg --name impactor-2025
   git remote add azure <git-url-from-above-command>
   git push azure main
   ```

4. **Monitor deployment:**
   ```bash
   az webapp log tail --resource-group impactor-rg --name impactor-2025
   ```

## Verify Deployment

- [ ] Check app logs show "Next.js is ready!"
- [ ] Check app logs show "Starting Streamlit..."
- [ ] Visit: `https://impactor-2025.azurewebsites.net`
- [ ] Verify background shader is visible
- [ ] Test all features work

## Files Required in Git

✅ **Must commit:**
- `start.sh` (startup script)
- `requirements.txt` (Python deps)
- `main.py` (Streamlit app)
- `src/` (all source code)
- `my-app/package.json` (Node.js deps)
- `my-app/next.config.ts` (Next.js config)
- `my-app/app/` (Next.js app code)
- `my-app/components/` (Next.js components)

❌ **Don't commit (in .gitignore):**
- `venv/` (Python virtual env)
- `my-app/node_modules/` (Node.js packages)
- `my-app/.next/` (Next.js build)
- `my-app/embed/` (Next.js static export)
- `__pycache__/` (Python cache)
- `.env` (secrets)

## Quick Commands

**View live logs:**
```bash
az webapp log tail --resource-group impactor-rg --name impactor-2025
```

**Restart app:**
```bash
az webapp restart --resource-group impactor-rg --name impactor-2025
```

**Stop app (save costs):**
```bash
az webapp stop --resource-group impactor-rg --name impactor-2025
```

**Start app:**
```bash
az webapp start --resource-group impactor-rg --name impactor-2025
```

**Delete everything:**
```bash
az group delete --name impactor-rg --yes
```

## Common Issues

### "Next.js not ready"
- Check Node.js version is 18+
- Verify `npm install --legacy-peer-deps` succeeds
- Check for port 3000 conflicts

### "Streamlit won't start"
- Verify `PORT=8000` is set
- Check `start.sh` has correct Streamlit flags
- Ensure Python 3.9+ runtime

### "Background not showing"
- Confirm both services running (check logs)
- Verify iframe src is `http://localhost:3000`
- Test Next.js URL directly

### "Deployment takes forever"
- First deploy installs all packages (5-10 min normal)
- Subsequent deploys are faster (~2-3 min)
- Check logs for any npm/pip errors

## App URL

After deployment: **https://impactor-2025.azurewebsites.net**

(Replace `impactor-2025` with your chosen app name)
