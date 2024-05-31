if [ ! -d ".venv" ]; then
    echo "setting up virtual environment and activating it"
    python3 -m venv .venv && source .venv/bin/activate
else
    echo "virtual environment already set, activating it"
    source .venv/bin/activate
fi
