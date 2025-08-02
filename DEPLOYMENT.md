# 🚀 Deploy Your Peptide Tag Generator

## Quick Deploy Options

### 1. Railway (Recommended - Free Tier)
**Best for:** Easy deployment with free tier

1. **Sign up** at [railway.app](https://railway.app)
2. **Connect your GitHub** repository
3. **Deploy** - Railway will automatically detect it's a Python app
4. **Get your URL** - Your app will be live at `https://your-app-name.railway.app`

### 2. Render (Popular - Free Tier)
**Best for:** Free hosting with good performance

1. **Sign up** at [render.com](https://render.com)
2. **Create a new Web Service**
3. **Connect your GitHub** repository
4. **Configure:**
   - **Build Command:** `pip install -r requirements.txt`
   - **Start Command:** `python app.py`
5. **Deploy** - Your app will be live at `https://your-app-name.onrender.com`

### 3. Heroku (Paid - $7/month minimum)
**Best for:** Established platform with good features

1. **Sign up** at [heroku.com](https://heroku.com)
2. **Install Heroku CLI**
3. **Deploy:**
   ```bash
   heroku create your-app-name
   git add .
   git commit -m "Initial deployment"
   git push heroku main
   ```

### 4. PythonAnywhere (Python-Specific)
**Best for:** Python-focused hosting

1. **Sign up** at [pythonanywhere.com](https://pythonanywhere.com)
2. **Upload your files** via Files tab
3. **Create a new Web app** (Flask)
4. **Configure WSGI file** to point to your app
5. **Deploy** - Your app will be live at `https://your-username.pythonanywhere.com`

## 🛠️ Local Testing Before Deployment

Test your app locally to make sure everything works:

```bash
# Install dependencies
pip install -r requirements.txt

# Run the app
python app.py

# Visit http://localhost:5000
```

## 📁 Files Included for Deployment

- ✅ `app.py` - Main Flask application
- ✅ `requirements.txt` - Python dependencies
- ✅ `Procfile` - Tells deployment platforms how to run the app
- ✅ `runtime.txt` - Specifies Python version
- ✅ `templates/` - HTML templates
- ✅ `static/` - CSS and static files

## 🔧 Environment Variables (Optional)

Some platforms let you set environment variables:

- `PORT` - Port number (usually set automatically)
- `FLASK_ENV` - Set to `production` for production

## 🌐 Custom Domain (Optional)

After deployment, you can add a custom domain:
1. **Railway/Render:** Settings → Custom Domains
2. **Heroku:** Settings → Domains
3. **PythonAnywhere:** Web tab → Add a new domain

## 📊 Monitoring & Analytics

Consider adding:
- **Google Analytics** for visitor tracking
- **Uptime monitoring** (UptimeRobot, Pingdom)
- **Error tracking** (Sentry)

## 🔒 Security Considerations

- ✅ HTTPS is automatically provided by most platforms
- ✅ Environment variables for sensitive data
- ✅ Input validation already implemented
- ✅ Rate limiting (consider adding for high traffic)

## 💰 Cost Comparison

| Platform | Free Tier | Paid Plans | Best For |
|----------|-----------|------------|----------|
| Railway | ✅ Yes | $5/month | Quick deployment |
| Render | ✅ Yes | $7/month | Good performance |
| Heroku | ❌ No | $7/month | Established platform |
| PythonAnywhere | ✅ Yes | $5/month | Python-specific |

## 🚀 Recommended: Railway

**Why Railway?**
- ✅ Free tier available
- ✅ Automatic deployments from GitHub
- ✅ Custom domains
- ✅ SSL certificates included
- ✅ Easy scaling
- ✅ Great developer experience

**Deploy in 5 minutes:**
1. Go to [railway.app](https://railway.app)
2. Sign up with GitHub
3. Click "New Project" → "Deploy from GitHub repo"
4. Select your repository
5. Deploy! 🎉

Your app will be live and accessible to everyone on the web! 