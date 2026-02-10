#!/bin/bash



# Check if Python3 and required packages are available
echo "📋 Checking Python environment..."
if ! command -v python3 &> /dev/null; then
    echo "❌ Python3 not found. Installing..."
    sudo apt update
    sudo apt install -y python3 python3-pip
fi

# Install required Python packages
echo "📦 Installing Python packages..."
pip3 install --break-system-packages flask flask-cors requests 2>/dev/null || \
pip3 install --user flask flask-cors requests 2>/dev/null || \
sudo apt install -y python3-flask python3-requests

# Create a standalone web interface
echo "🔧 Creating standalone web interface..."
cat > antibody_web.py << 'EOF'
#!/usr/bin/env python3
import re
import json
import requests
import subprocess
from flask import Flask, render_template_string, request, jsonify
from flask_cors import CORS
import os
import sys

app = Flask(__name__)
CORS(app)

MODEL_NAME = "hf.co/silicobio/peleke-phi-4-gguf:Q4_K_M"

# HTML template embedded in Python
HTML_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Antibody Designer</title>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { 
            font-family: Arial, sans-serif; 
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px; 
            min-height: 100vh;
        }
        .container { max-width: 1200px; margin: 0 auto; }
        .header { 
            background: white; 
            padding: 30px; 
            border-radius: 15px; 
            text-align: center; 
            margin-bottom: 30px; 
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }
        .header h1 { color: #2c3e50; font-size: 2.5rem; margin-bottom: 10px; }
        .header p { color: #7f8c8d; font-size: 1.2rem; }
        .main-content { display: grid; grid-template-columns: 1fr 1fr; gap: 30px; }
        .card { 
            background: white; 
            padding: 25px; 
            border-radius: 15px; 
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
            transition: transform 0.2s;
        }
        .card:hover { transform: translateY(-5px); }
        .card h2 { color: #2c3e50; margin-bottom: 20px; }
        .form-group { margin-bottom: 20px; }
        .form-label { display: block; margin-bottom: 8px; font-weight: bold; color: #34495e; }
        .form-input { 
            width: 100%; 
            padding: 12px; 
            border: 2px solid #ecf0f1; 
            border-radius: 8px; 
            font-size: 1rem; 
            min-height: 120px;
            font-family: 'Courier New', monospace;
        }
        .form-input:focus { outline: none; border-color: #3498db; }
        .btn { 
            background: linear-gradient(135deg, #3498db, #2980b9); 
            color: white; 
            border: none; 
            padding: 12px 24px; 
            border-radius: 8px; 
            cursor: pointer; 
            font-size: 1rem; 
            font-weight: bold;
            transition: all 0.3s;
        }
        .btn:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(52, 152, 219, 0.4); }
        .btn:disabled { opacity: 0.6; cursor: not-allowed; transform: none; }
        .loading { display: none; margin-top: 15px; color: #7f8c8d; font-style: italic; }
        .result-card { 
            background: #ecf0f1; 
            border: 2px solid #27ae60; 
            border-radius: 10px; 
            padding: 20px; 
            margin-top: 20px; 
            display: none; 
        }
        .chain-section { margin-bottom: 15px; }
        .chain-title { font-weight: bold; color: #2c3e50; margin-bottom: 8px; }
        .chain-sequence { 
            background: white; 
            padding: 10px; 
            border-radius: 5px; 
            font-family: 'Courier New', monospace; 
            font-size: 0.9rem; 
            word-break: break-all;
            border-left: 4px solid #3498db;
        }
        .error { 
            background: #e74c3c; 
            color: white; 
            padding: 15px; 
            border-radius: 8px; 
            margin-top: 15px; 
            display: none; 
        }
        .example { 
            background: #d5f4e6; 
            border: 1px solid #27ae60; 
            padding: 15px; 
            border-radius: 8px; 
            margin-top: 15px; 
        }
        .status { 
            position: fixed; 
            top: 20px; 
            right: 20px; 
            padding: 10px 15px; 
            border-radius: 20px; 
            font-weight: bold;
        }
        .status-healthy { background: #d5f4e6; color: #27ae60; }
        .status-unhealthy { background: #fadbd8; color: #e74c3c; }
        @media (max-width: 768px) {
            .main-content { grid-template-columns: 1fr; }
            .header h1 { font-size: 2rem; }
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Peleke Antibody Designer</h1>
            <p>Generate antibody sequences from antigen epitopes using AI</p>
        </div>

        <div id="status" class="status status-healthy">✅ System Ready</div>

        <div class="main-content">
            <div class="card">
                <h2>🧬 Input Sequence</h2>
                
                <div class="form-group">
                    <label class="form-label">Antigen Sequence (mark epitope regions with brackets)</label>
                    <textarea id="antigenSequence" class="form-input" 
                              placeholder="Enter your antigen sequence with epitope regions marked as [A], [B], etc.
Example: NPPTFSPALL[S][F][V]LNWY"></textarea>
                </div>

                <button id="generateBtn" class="btn" onclick="generateAntibody()">
                    🧪 Generate Antibody
                </button>

                <div id="loading" class="loading">⏳ Generating antibody sequence...</div>

                <div class="example">
                    <strong>Example Input:</strong><br>
                    <code>NPPTFSPALL[S][F][V]LNWY</code>
                </div>

                <div id="error" class="error"></div>
            </div>

            <div class="card">
                <h2>⚡ Generated Antibody</h2>

                <div id="results" class="result-card">
                    <div class="chain-section">
                        <div class="chain-title">🔗 Heavy Chain</div>
                        <div id="heavyChain" class="chain-sequence"></div>
                    </div>
                    <div class="chain-section">
                        <div class="chain-title">🔗 Light Chain</div>
                        <div id="lightChain" class="chain-sequence"></div>
                    </div>
                    <div class="chain-section">
                        <div class="chain-title">📋 Full Sequence</div>
                        <div id="fullSequence" class="chain-sequence"></div>
                    </div>
                </div>

                <div style="margin-top: 20px; color: #7f8c8d; font-style: italic;">
                    💡 Results will appear here after generating an antibody sequence
                </div>
            </div>
        </div>
    </div>

    <script>
        async function generateAntibody() {
            const antigenSequence = document.getElementById('antigenSequence').value.trim();
            const errorDiv = document.getElementById('error');
            const resultsDiv = document.getElementById('results');
            const loadingDiv = document.getElementById('loading');
            const generateBtn = document.getElementById('generateBtn');

            errorDiv.style.display = 'none';
            resultsDiv.style.display = 'none';

            if (!antigenSequence) {
                showError('Please enter an antigen sequence');
                return;
            }

            if (!/\\[[A-Z]\\]/.test(antigenSequence)) {
                showError('Please mark epitope regions with brackets, e.g., [S][F][V]');
                return;
            }

            loadingDiv.style.display = 'block';
            generateBtn.disabled = true;

            try {
                const response = await fetch('/api/generate', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ antigen_sequence: antigenSequence })
                });

                const data = await response.json();

                if (!response.ok) {
                    throw new Error(data.error || 'Generation failed');
                }

                document.getElementById('heavyChain').textContent = data.heavy_chain || 'Not available';
                document.getElementById('lightChain').textContent = data.light_chain || 'Not available';
                document.getElementById('fullSequence').textContent = data.antibody_sequence || 'No sequence generated';

                resultsDiv.style.display = 'block';

            } catch (error) {
                showError(`Error: ${error.message}`);
            } finally {
                loadingDiv.style.display = 'none';
                generateBtn.disabled = false;
            }
        }

        function showError(message) {
            const errorDiv = document.getElementById('error');
            errorDiv.textContent = message;
            errorDiv.style.display = 'block';
        }
    </script>
</body>
</html>
'''

def format_prompt(antigen_sequence):
    """Convert bracket notation [A] to <epi>A</epi> format"""
    epitope_seq = re.sub(r'\[([A-Z])\]', r'<epi>\1</epi>', antigen_sequence)
    formatted_str = f"Antigen: {epitope_seq}<|im_end|>\nAntibody:"
    return formatted_str

def generate_antibody(antigen_sequence):
    """Generate antibody sequence using Ollama"""
    try:
        prompt = format_prompt(antigen_sequence)
        
        # Try direct ollama command first
        try:
            result = subprocess.run([
                'ollama', 'run', MODEL_NAME, prompt
            ], capture_output=True, text=True, timeout=120)
            
            if result.returncode == 0:
                return result.stdout.strip().replace('Antibody:', '').strip()
        except:
            pass
        
        # Fall back to API call
        payload = {
            "model": MODEL_NAME,
            "prompt": prompt,
            "stream": False,
            "options": {
                "temperature": 0.7,
                "num_predict": 1000,
                "stop": ["<|im_end|>", "\n\n"]
            }
        }
        
        response = requests.post("http://localhost:11434/api/generate", json=payload, timeout=120)
        response.raise_for_status()
        
        result = response.json()
        antibody_response = result.get("response", "").strip()
        return antibody_response.replace("Antibody:", "").strip()
        
    except Exception as e:
        return f"Error: {str(e)}"

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/api/generate', methods=['POST'])
def api_generate():
    try:
        data = request.get_json()
        antigen_sequence = data.get('antigen_sequence', '').strip()
        
        if not antigen_sequence:
            return jsonify({'error': 'Antigen sequence is required'}), 400
        
        if not re.search(r'\[[A-Z]\]', antigen_sequence):
            return jsonify({'error': 'Please mark epitope regions with brackets, e.g., [S][F][V]'}), 400
        
        antibody_sequence = generate_antibody(antigen_sequence)
        
        chains = antibody_sequence.split('|') if '|' in antibody_sequence else [antibody_sequence]
        heavy_chain = chains[0].strip() if len(chains) > 0 else ""
        light_chain = chains[1].strip() if len(chains) > 1 else ""
        
        return jsonify({
            'antigen_sequence': antigen_sequence,
            'antibody_sequence': antibody_sequence,
            'heavy_chain': heavy_chain,
            'light_chain': light_chain
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/health')
def health_check():
    try:
        # Check if ollama is available
        result = subprocess.run(['ollama', 'list'], capture_output=True)
        if result.returncode == 0:
            return jsonify({'status': 'healthy', 'ollama_running': True})
        else:
            return jsonify({'status': 'unhealthy', 'ollama_running': False}), 503
    except:
        return jsonify({'status': 'unhealthy', 'ollama_running': False}), 503

if __name__ == '__main__':
    print("🧬 Starting Antibody Designer Web Interface...")
    print("📝 Access at: http://localhost:8081")
    print("🔧 Press Ctrl+C to stop")
    app.run(host='0.0.0.0', port=8081, debug=False)
EOF

# Make it executable
chmod +x antibody_web.py

# Check if Ollama is running and has the model
echo "🔍 Checking Ollama status..."
if ! command -v ollama &> /dev/null; then
    echo "❌ Ollama not found in PATH"
    echo "   Make sure Ollama is installed and accessible"
    exit 1
fi

if ! ollama list | grep -q "hf.co/silicobio/peleke-phi-4-gguf:Q4_K_M"; then
    echo "📥 Model not found. Downloading..."
    ollama pull hf.co/silicobio/peleke-phi-4-gguf:Q4_K_M
fi

echo "✅ Setup complete!"
echo ""
echo "🚀 Starting Antibody Designer..."
echo "📝 Access at: http://localhost:8081"
echo "🔧 Press Ctrl+C to stop"
echo ""

# Start the web interface
python3 antibody_web.py
