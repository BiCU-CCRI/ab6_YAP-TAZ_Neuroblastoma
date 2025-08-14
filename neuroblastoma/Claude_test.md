Use https://github.com/imnmv/clauder

## Install nodejs first.
Run this in the terminal.
- Download and install nvm:
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.3/install.sh | bash
- In lieu of restarting the shell `\. "$HOME/.nvm/nvm.sh"`
- Download and install Node.js: `nvm install 22`
- Install claude `npm install -g @anthropic-ai/claude-code`
- In R install and activate claude: ```devtools::install_github("IMNMV/ClaudeR")
library(ClaudeR)
install_cli(tools = "claude")```
